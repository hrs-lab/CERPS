"""3D space: cell positions, mechanical relaxation, contact inhibition
and a paracrine signalling field.

* Cells occupy spheres of radius ``cell_radius``; overlapping pairs are
  relaxed apart by a fraction of the overlap per step (overdamped
  soft-sphere mechanics; Van Liedekerke et al. 2018; Osborne et al. 2017).
* Divisions place the two daughters symmetrically about the mother along a
  random axis.
* Local density (neighbour count within the inhibition radius) feeds back
  on the proliferation rate (contact inhibition; Shraiman's mechanical
  growth control; Zimmermann et al. 2016).
* Optional paracrine field c(r, t): diffusion-decay on a regular grid
  with cell-based secretion (PhysiCell/PhysiBoSS-style; Letort et al.
  2019); cells sense the local field, which adds to gene activation.
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
        """One step of overdamped overlap relaxation + boundary."""
        from scipy.spatial import cKDTree
        pos = self.pos
        if len(pos) < 2:
            return pos
        tree = cKDTree(pos)
        pairs = tree.query_pairs(2.0 * self.R, output_type="ndarray")
        if len(pairs):
            d = pos[pairs[:, 1]] - pos[pairs[:, 0]]
            dist = np.linalg.norm(d, axis=1)
            dist = np.maximum(dist, 1e-9)
            overlap = 2.0 * self.R - dist
            m = overlap > 0
            if m.any():
                push = (0.5 * self.cfg.relax * overlap[m])[:, None] * \
                    d[m] / dist[m, None]
                # symmetric displacement
                np.add.at(pos, pairs[m, 0], -push)
                np.add.at(pos, pairs[m, 1], push)
        # motility
        if self.cfg.motility > 0:
            raise NotImplementedError("motility handled by simulator")
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
        """One explicit diffusion-decay step + cell-based sources.

        secretion: (n_cells,) source strengths (e.g. stemness-dependent).
        Updates self.field (3D grid); returns (n_cells,) sampled field.
        """
        cfg = self.cfg
        if self.field is None:
            self.field = np.zeros(self._grid_shape)
        f = self.field
        h = cfg.grid_h
        # sources (nearest grid point; grid is fine enough for demos)
        idx = self._to_grid(self.pos)
        np.add.at(f, (idx[:, 0], idx[:, 1], idx[:, 2]), secretion * dt)
        # decay
        f *= (1.0 - cfg.decay_rate * dt)
        # diffusion (FTCS, von Neumann stable for D dt/h^2 <= 1/6)
        c = cfg.diff_rate * dt / h ** 2
        c = min(c, 1.0 / 6.0)
        lap = -6.0 * f
        for ax in range(3):
            lap += np.roll(f, 1, ax) + np.roll(f, -1, ax)
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
