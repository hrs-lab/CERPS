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
        if cfg.paracrine:
            h = cfg.grid_h
            self._grid_shape = (int(2 * self.R_dom / h) + 1,) * 3
            self._grid_origin = -self.R_dom
            self.field = None

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
        idx = self._to_grid(self.pos)
        for _ in range(n_sub):
            # sources (nearest grid point)
            np.add.at(f, (idx[:, 0], idx[:, 1], idx[:, 2]),
                      secretion * sdt * src_scale)
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
        self.field = f
        return self.sample(self.pos)

    field = None

    def _to_grid(self, pos):
        h = self.cfg.grid_h
        idx = ((pos - self._grid_origin) / h).astype(int)
        shape = np.array(self._grid_shape)
        return np.clip(idx, 0, shape - 1)

    def sample(self, pos):
        """Nearest-neighbour sampling of the field at cell positions."""
        if self.field is None:
            return np.zeros(len(pos))
        idx = self._to_grid(pos)
        return self.field[idx[:, 0], idx[:, 1], idx[:, 2]]
