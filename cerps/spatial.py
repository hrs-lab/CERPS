"""3D space: cell positions, mechanical relaxation, contact inhibition
and a paracrine signalling field.

* Cells occupy spheres of radius ``cell_radius``; overlapping pairs are
  relaxed apart with an exponentially decaying overlap (overdamped
  soft-sphere mechanics; Van Liedekerke et al. 2018/2019/2020).
* Divisions place the two daughters symmetrically about the mother along a
  random axis.
* Local density (neighbour count within the inhibition radius) feeds back
  on the proliferation rate (contact inhibition; Shraiman's mechanical
  growth control; Zimmermann et al. 2016).
* Optional paracrine field c(r, t): diffusion-decay on a regular grid
  with cell-based secretion (PhysiCell/PhysiBoSS-style; Letort et al.
  2019; Du, Nie & Holmes 2015); cells sense the local field, which adds
  to gene activation.  The solver sub-cycles to respect the 3-D CFL limit
  Fo = D*dt/h^2 <= 1/6 (never clamping the nominal diffusivity) and uses
  reflecting boundaries consistent with the spherical cell domain.
"""
import numpy as np


class Space:
    def __init__(self, cfg, domain=None):
        self.cfg = cfg
        self.R = cfg.cell_radius
        self.R_dom = cfg.domain_radius
        self.pos = None            # (n, 3)
        self._grid = None
        self._grid_shape = None
        self._grid_origin = None
        self._sphere_mask = None
        if cfg.paracrine:
            h = cfg.grid_h
            self._grid_shape = (int(2 * self.R_dom / h) + 1,) * 3
            self._grid_origin = -self.R_dom
            self.field = None
            if getattr(cfg, "sphere_mask", False):
                ax = self._grid_origin + h * np.arange(self._grid_shape[0])
                X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
                self._sphere_mask = (X ** 2 + Y ** 2 + Z ** 2) <= \
                    self.R_dom ** 2

    # ------------------------------------------------------------------
    def init_positions(self, n, rng, cluster_radius=None):
        """Random positions inside a small cluster (default: sphere of
        radius ~ n^(1/3) * R)."""
        if cluster_radius is None:
            cluster_radius = max(2.0 * self.R,
                                 1.2 * self.R * n ** (1.0 / 3.0))
        pos = rng.normal(size=(n, 3))
        pos *= cluster_radius / np.linalg.norm(pos, axis=1,
                                               keepdims=True) * \
            rng.random(n)[:, None] ** (1.0 / 3.0)
        self.pos = pos
        return pos

    def add(self, new_pos):
        self.pos = np.vstack([self.pos, new_pos])

    def remove(self, keep_mask):
        self.pos = self.pos[keep_mask]

    # ------------------------------------------------------------------
    def local_density(self, radius_mult=None):
        """Neighbour counts within the inhibition radius (excluding self)."""
        from scipy.spatial import cKDTree
        rm = radius_mult if radius_mult is not None else \
            self.cfg.inhibition_radius
        tree = cKDTree(self.pos)
        counts = np.array(tree.query_ball_point(
            self.pos, rm * 2 * self.R, return_length=True))
        return np.maximum(counts - 1, 0)

    def relax(self, dt):
        """One step of overdamped overlap relaxation + boundary.

        The overlap decays exponentially with the mechanical relaxation
        time ``relax_time`` [h]:  overlap <- overlap * exp(-dt/tau).
        This makes the mechanics time-step consistent (the previous
        version applied a fixed fraction per *step*, so changing dt
        changed the mechanical rate without changing the dynamics it
        converged to).

        ``relax`` (fraction per step at a nominal dt of 0.25 h) is kept
        for backward compatibility: if ``relax_time`` is None it is
        converted as tau = -0.25 / ln(1 - relax), which reproduces the
        old behaviour exactly at dt = 0.25 h.
        """
        from scipy.spatial import cKDTree
        pos = self.pos
        if len(pos) < 2:
            return pos
        tau = self.cfg.relax_time
        if tau is None:
            r = min(max(self.cfg.relax, 1e-6), 0.999999)
            tau = -0.25 / np.log(1.0 - r)
        frac = 1.0 - np.exp(-dt / max(tau, 1e-9))
        tree = cKDTree(pos)
        pairs = tree.query_pairs(2.0 * self.R, output_type="ndarray")
        if len(pairs):
            d = pos[pairs[:, 1]] - pos[pairs[:, 0]]
            dist = np.linalg.norm(d, axis=1)
            dist = np.maximum(dist, 1e-9)
            overlap = 2.0 * self.R - dist
            m = overlap > 0
            if m.any():
                push = (0.5 * frac * overlap[m])[:, None] * \
                    d[m] / dist[m, None]
                # symmetric displacement
                np.add.at(pos, pairs[m, 0], -push)
                np.add.at(pos, pairs[m, 1], push)
        # motility is handled by the simulator
        # reflecting spherical boundary
        if self.R_dom is not None and self.R_dom > 0:
            r = np.linalg.norm(pos, axis=1)
            out = r > self.R_dom - self.R
            if out.any():
                pos[out] *= ((self.R_dom - self.R) /
                             r[out])[:, None]
        return pos

    # ------------------------------------------------------------------
    def division_positions(self, idx, rng):
        """Positions for the two daughters of cells at ``idx``."""
        axis = rng.normal(size=(len(idx), 3))
        axis /= np.maximum(np.linalg.norm(axis, axis=1, keepdims=True),
                           1e-9)
        d = self.R * 0.6
        p1 = self.pos[idx] - axis * d
        p2 = self.pos[idx] + axis * d
        return p1, p2

    # ------------------------------------------------------------------
    # paracrine field
    # ------------------------------------------------------------------
    def field_step(self, secretion, dt):
        """Advance the paracrine field by ``dt`` and sample it at the cells.

        Solves  dc/dt = D lap(c) - k c + sources  on the regular grid with
        an explicit FTCS scheme.  Two properties matter for quantitative
        use:

        1. **No silent CFL clamping.**  The 3-D von Neumann stability limit
           is Fo = D*dt/h^2 <= 1/6.  Instead of clamping the diffusion
           number (which silently replaces the nominal diffusivity by a
           smaller effective one -- e.g. D=25, dt=0.25, h=4 would run at
           D_eff = 10.67), the step is sub-cycled so that the *nominal* D
           is honoured and the scheme stays stable.
        2. **Reflecting (zero-flux) boundaries**, matching the reflecting
           spherical cell domain.  The previous implementation used
           ``np.roll``, i.e. periodic boundaries, which leaked ~8% of the
           signal across the domain.

        ``secretion`` : (n_cells,) source strengths (e.g. stemness-
            dependent).  Units: the field carries secretion*time in
            arbitrary units; ``SpaceConfig.response`` maps it onto gene
            activation.  Set ``normalize_source=True`` to divide the
            injection by h^3, giving a proper concentration field.

        Returns (n_cells,) sampled field.
        """
        cfg = self.cfg
        if self.field is None:
            self.field = np.zeros(self._grid_shape)
        f = self.field
        h = cfg.grid_h
        D, k = cfg.diff_rate, cfg.decay_rate
        # --- adaptive sub-stepping (no CFL clamp) -----------------------
        cfl = D * dt / h ** 2
        n_sub = max(1, int(np.ceil(cfl / (1.0 / 6.0) - 1e-12)))
        sdt = dt / n_sub
        c = D * sdt / h ** 2                      # <= 1/6 by construction
        src_scale = (1.0 / h ** 3) if cfg.normalize_source else 1.0
        linear = getattr(cfg, "field_interp", "nearest") == "linear"
        idx = self._to_grid(self.pos)
        if linear:
            # precompute the 8 trilinear corner indices/weights (position is
            # fixed during the sub-steps)
            ci0, ci1, cw0, cw1 = self._tri(self.pos)
        for _ in range(n_sub):
            q = secretion * sdt * src_scale
            if linear:
                # conservative trilinear deposition: the 8 weights sum to 1,
                # so total injected mass is unchanged, but the source is no
                # longer a grid-aligned point (removes O(h) anisotropy)
                for aa in (0, 1):
                    ia = ci1[:, 0] if aa else ci0[:, 0]
                    wa = cw1[:, 0] if aa else cw0[:, 0]
                    for bb in (0, 1):
                        ib = ci1[:, 1] if bb else ci0[:, 1]
                        wb = cw1[:, 1] if bb else cw0[:, 1]
                        for cc in (0, 1):
                            ic = ci1[:, 2] if cc else ci0[:, 2]
                            wc = cw1[:, 2] if cc else cw0[:, 2]
                            np.add.at(f, (ia, ib, ic), q * wa * wb * wc)
            else:
                # sources (nearest grid point)
                np.add.at(f, (idx[:, 0], idx[:, 1], idx[:, 2]), q)
            # decay (explicit Euler, stable for k*sdt < 2)
            f *= (1.0 - min(k * sdt, 1.0))
            # diffusion, FTCS with reflecting boundaries (edge padding)
            fp = np.pad(f, 1, mode="edge")
            # explicit 6-neighbour sum from the padded array
            lap = (-6.0 * f
                   + fp[2:, 1:-1, 1:-1] + fp[:-2, 1:-1, 1:-1]
                   + fp[1:-1, 2:, 1:-1] + fp[1:-1, :-2, 1:-1]
                   + fp[1:-1, 1:-1, 2:] + fp[1:-1, 1:-1, :-2])
            f += c * lap
            if self._sphere_mask is not None:
                # Dirichlet boundary at r = R_dom: no field outside the cell
                # domain (removes the cell-free cube-corner reservoir)
                f[~self._sphere_mask] = 0.0
        self.field = f
        return self.sample(self.pos)

    field = None

    def _to_grid(self, pos):
        h = self.cfg.grid_h
        idx = ((pos - self._grid_origin) / h).astype(int)
        shape = np.array(self._grid_shape)
        return np.clip(idx, 0, shape - 1)

    def _tri(self, pos):
        """Trilinear cell indices and weights for ``pos`` (cube grid).

        Returns ``(i0, i1, w0, w1)``, each of shape (m, 3): the lower and
        upper grid index along every axis (clipped to the grid) and the
        corresponding weights (``w0 + w1 == 1``).  At the upper boundary
        ``i0 == i1`` and ``w1 == 0``, so the 8-point stencil degenerates
        gracefully to the face node.
        """
        h = self.cfg.grid_h
        n = np.array(self._grid_shape)
        f = np.clip((pos - self._grid_origin) / h, 0.0, n[None, :] - 1.0)
        i0 = np.floor(f).astype(int)
        i1 = np.minimum(i0 + 1, n[None, :] - 1)
        w1 = f - i0
        return i0, i1, 1.0 - w1, w1

    def sample(self, pos):
        """Sample the field at cell positions.

        ``field_interp="nearest"`` (default) reproduces the validated
        baseline (piecewise-constant, O(h) error).  ``"linear"`` uses
        trilinear interpolation, which is exact for a field that is linear
        over one grid cell and removes the grid anisotropy.
        """
        if self.field is None:
            return np.zeros(len(pos))
        if getattr(self.cfg, "field_interp", "nearest") == "linear":
            i0, i1, w0, w1 = self._tri(pos)
            f = self.field
            out = np.zeros(len(pos))
            for aa in (0, 1):
                ia = i1[:, 0] if aa else i0[:, 0]
                wa = w1[:, 0] if aa else w0[:, 0]
                for bb in (0, 1):
                    ib = i1[:, 1] if bb else i0[:, 1]
                    wb = w1[:, 1] if bb else w0[:, 1]
                    for cc in (0, 1):
                        ic = i1[:, 2] if cc else i0[:, 2]
                        wc = w1[:, 2] if cc else w0[:, 2]
                        out += wa * wb * wc * f[ia, ib, ic]
            return out
        idx = self._to_grid(pos)
        return self.field[idx[:, 0], idx[:, 1], idx[:, 2]]
