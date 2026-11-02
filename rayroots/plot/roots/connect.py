# Frozen method (v1 + v2). Do not edit here: a change is a new version.
"""v2 (root connectivity): own fragments with measured thickness, and multi-stem growth with structural rules.
Fragments ("block fragments"): points voxelised at VOX, space cut into cubes of BLOCK voxels (L = 0.20 m); a
fragment is a 26-connected piece of occupied voxels inside one cube. Fragments are independent of any seed (unlike
RCT's geodesic level sets, which are built from the stems through the mat). Fragments touching across a cube face
(26-neighbour voxels) are graph neighbours. Pieces with < MIN_PTS points are merged into the neighbour they touch most.
Per fragment from its own points: centroid, PCA axis, length (p2-p98 along the axis), width (p5-p95 along the second
axis; for any shape this is <= the real diameter, so it is a lower bound: a sheet or a mat lump gets a large width,
a root stub cut across by a cube corner gets a small one), and quality: 'measured' when n >= N_MEAS and
length >= ELONG x width (axis is the real axis, so width ~ diameter of the visible part).
Growth: see Growth; every stop is logged with its rule. thick_frags(): fragments for the second pass."""
import heapq
from collections import defaultdict

import numpy as np
from scipy import sparse
from scipy.sparse import csgraph

VOX, BLOCK = 0.025, 8          # 2.5 cm voxels, 0.20 m cubes: >= 2x a prop root's diameter (axis measurable),
                               # short enough that an arching root (radius >= 0.5 m) bends < 25 deg inside one
MIN_PTS = 10                   # smaller pieces are merged into their best-touching neighbour
N_MEAS, ELONG = 30, 1.5        # measured thickness: enough points and clearly elongated
# Growth rules (prototype values, not tuned)
TOL, ABS_TOL = 0.30, 0.01      # child width <= bound x (1 + TOL) + ABS_TOL (occlusion changes the visible chord)
BOUND_K = 3                    # bound = largest of the last BOUND_K measured widths on the path (robust to one
                               # under-measured fragment)
SPLIT_TOL = 0.30               # children: sum of measured width^2 <= bound^2 x (1 + SPLIT_TOL)^2 (pipe model)
TURN_MAX = 40.0                # deg between axes of consecutive measured pieces (~0.2 m): curvature radius >= ~0.3 m
STEP_MAX = 60.0                # deg: the centroid step to any child goes forward (cube pieces' centroids are noisy)
FAN_N, FAN_R = 4, 0.6          # fan (mat): >= FAN_N candidate neighbours whose directions spread (resultant < FAN_R)
LOOK = 2                       # coherence: in a fan (and at the stem -> root step) a child must carry on LOOK steps
BRANCH_MAX = 90.0              # deg: at a non-fan split a measured, coherent child may leave at up to this angle
                               # (branching is not bending; biological, fixed, not tuned; Jasper 2026-10-18).
                               # None = strict (prototype)
FINE = 4                       # second pass: fragments that stopped growth on thickness are re-cut in cubes of
                               # FINE voxels (0.10 m), which separates roots touching inside one 0.20 m cube
R_REACH = 7.0                  # m: horizontal reach from the stem
TIE_REL, TIE_ABS = 0.10, 0.10  # near-tie: second tree's valid cost <= best + max(TIE_REL x best, TIE_ABS)
OFFS = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)
                 if (a, b, c) > (0, 0, 0)])           # 13 forward neighbour offsets


# %% Fragments
def fragments(xyz, fine=None):
    """Fragment id per point, n fragments, adjacent pairs (a < b) and their contacts. fine: bool per point; voxels
    holding such points are cut in cubes of FINE voxels instead of BLOCK (second pass)."""
    lo = xyz.min(0) - VOX
    ijk = np.floor((xyz - lo) / VOX).astype(np.int64)
    dims = ijk.max(0) + 2
    key = (ijk[:, 0] * dims[1] + ijk[:, 1]) * dims[2] + ijk[:, 2]
    vk, inv = np.unique(key, return_inverse=True)
    vijk = np.c_[vk // (dims[1] * dims[2]), (vk // dims[2]) % dims[1], vk % dims[2]]
    nv = len(vk)
    vfine = np.zeros(nv, bool) if fine is None else np.bincount(inv, fine, nv) > 0
    cube = vijk // np.where(vfine, FINE, BLOCK)[:, None]
    intra, inter = [], []
    for o in OFFS:
        q = vijk + o
        ok = (q >= 0).all(1) & (q < dims).all(1)
        nk = (q[:, 0] * dims[1] + q[:, 1]) * dims[2] + q[:, 2]
        pos = np.minimum(np.searchsorted(vk, nk), nv - 1)
        hit = ok & (vk[pos] == nk)
        a, b = np.flatnonzero(hit), pos[hit]
        same = (vfine[a] == vfine[b]) & (cube[a] == cube[b]).all(1)
        intra.append(np.c_[a[same], b[same]]); inter.append(np.c_[a[~same], b[~same]])
    e = np.vstack(intra)
    g = sparse.coo_matrix((np.ones(len(e), np.int8), (e[:, 0], e[:, 1])), shape=(nv, nv))
    nf, vfrag = csgraph.connected_components(g, directed=False)
    e = np.vstack(inter)
    fa, fb = vfrag[e[:, 0]], vfrag[e[:, 1]]
    frag = vfrag[inv]
    npts = np.bincount(frag, minlength=nf)
    # merge small pieces into the neighbour they touch most (one pass; a small piece without neighbours stays)
    pairs = np.c_[np.r_[fa, fb], np.r_[fb, fa]]
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    u, cnt = np.unique(pairs, axis=0, return_counts=True)
    small = npts[u[:, 0]] < MIN_PTS
    u, cnt = u[small], cnt[small]
    big_nb = npts[u[:, 1]] >= MIN_PTS
    order = np.lexsort((-cnt, ~big_nb, u[:, 0]))      # per small piece: big neighbours first, most contacts
    u = u[order]
    first = np.r_[True, u[1:, 0] != u[:-1, 0]]
    remap = np.arange(nf)
    remap[u[first, 0]] = u[first, 1]
    for _ in range(3):                                  # small -> small -> big chains
        remap = remap[remap]
    _, remap = np.unique(remap, return_inverse=True)
    frag = remap[frag]; fa, fb = remap[fa], remap[fb]
    nf = remap.max() + 1
    k = fa != fb
    a, b = np.minimum(fa[k], fb[k]), np.maximum(fa[k], fb[k])
    u, cnt = np.unique(np.c_[a, b], axis=0, return_counts=True)
    return frag, nf, u, cnt


def _group_quantile(vals, grp, starts, n, q):
    o = np.lexsort((vals, grp))
    idx = starts + np.round(q * (n - 1)).astype(np.int64)
    return vals[o][idx]


def measure(xyz, frag, nf):
    """Per fragment: n, centroid, axis (unit), length, width, depth (third axis), measured flag."""
    n = np.bincount(frag, minlength=nf).astype(np.float64)
    c = np.stack([np.bincount(frag, xyz[:, i], nf) for i in range(3)], 1) / n[:, None]
    d = xyz - c[frag]
    cov = np.zeros((nf, 3, 3))
    for i in range(3):
        for j in range(i, 3):
            cov[:, i, j] = cov[:, j, i] = np.bincount(frag, d[:, i] * d[:, j], nf) / n
    w, v = np.linalg.eigh(cov)                      # ascending
    e1, e2, e3 = v[:, :, 2], v[:, :, 1], v[:, :, 0]
    starts = np.r_[0, np.cumsum(n)[:-1]].astype(np.int64)
    ni = n.astype(np.int64)
    out = {"n": ni, "c": c, "axis": e1}
    for name, e, (qa, qb) in (("length", e1, (0.02, 0.98)), ("width", e2, (0.05, 0.95)), ("depth", e3, (0.05, 0.95))):
        p = np.einsum("nd,nd->n", d, e[frag])
        out[name] = _group_quantile(p, frag, starts, ni, qb) - _group_quantile(p, frag, starts, ni, qa)
    out["measured"] = (ni >= N_MEAS) & (out["length"] >= ELONG * out["width"])
    return out


def adjacency(nf, pairs):
    nb = [[] for _ in range(nf)]
    for a, b in pairs:
        nb[a].append(b); nb[b].append(a)
    return nb


# %% Growth
def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def _angle(a, b):
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))))


class Growth:
    """All stems at once, best-first on path cost. State per (tree, fragment): cost, parent, bound history (last
    measured widths), direction. A fragment is settled by the first tree that pops it (lowest cost); a later valid
    arrival of another tree within the tie margin marks it a near-tie.
    Direction: leaving a measured fragment = its axis (oriented along the path); leaving an unmeasured one = the
    incoming direction nudged by the centroid step (centroid steps between cube pieces are noisy).
    Rules when a settled fragment P of tree S is expanded, over its candidate neighbours (not S's own, not another
    tree's RCT fragments, within R_REACH of S's stem):
      thickness  lower bound of the child's diameter (width if measured, else its smallest extent: unknown never
                 counts as thicker) <= bound x (1 + TOL) + ABS_TOL; unmeasured children carry the bound on
      direction  measured child: its axis turns <= TURN_MAX from P's direction ('bend'); every child: the centroid
                 step lies within STEP_MAX of P's direction ('backward'); up or down alike
      fan        >= FAN_N short (unmeasured) pieces around P spreading in all directions: only the best-aligned
                 passing child that carries on coherently (LOOK steps) continues; the rest stop ('mat')
      split      passing children in order of turn; measured widths^2 summed <= bound^2 (1 + SPLIT_TOL)^2
    The stem -> first root step is exempt from thickness (bound = stem diameter) and direction, but a first root
    fragment must carry on coherently (LOOK steps)."""

    def __init__(self, f, nb, seeds, stem_diam, stem_xy, rct_owner):
        self.f, self.nb = f, nb
        self.seeds, self.stem_diam, self.stem_xy, self.rct_owner = seeds, stem_diam, stem_xy, rct_owner
        nf = len(f["n"])
        self.owner = np.full(nf, -1, np.int64)
        self.cost = np.full(nf, np.inf)
        self.parent = np.full(nf, -1, np.int64)
        self.bound = np.full(nf, np.nan)
        self.tie = np.zeros(nf, bool)
        self.tie_with = np.full(nf, -1, np.int64)
        self.is_seed = np.zeros(nf, bool)
        self.stop = {}                                # (tree, frag) -> (rule, from fragment, value, bound)
        self.late = {}                                # (tree, frag) -> lowest valid cost of a tree that came
                                                      # after the owner (lost the competition)
        self.heap = []
        self.lower = np.where(f["measured"], f["width"], f["depth"])

    def _ok_cand(self, s, c):
        ro = self.rct_owner[c]
        if ro >= 0 and ro != s:
            self.stop.setdefault((s, c), ("other tree (RCT)", -1, np.nan, np.nan))
            return False
        if np.linalg.norm(self.f["c"][c, :2] - self.stem_xy[s]) > R_REACH:
            self.stop.setdefault((s, c), ("reach", -1, np.nan, np.nan))
            return False
        return True

    def _out_dir(self, c, din, step):
        """Direction leaving c, reached with direction din (None at the stem) by centroid step."""
        ref = step if din is None else din
        if self.f["measured"][c]:
            a = self.f["axis"][c]
            return a if np.dot(a, ref) >= 0 else -a
        return step if din is None else _unit(2 * din + step)

    def _turn(self, c, din, step):
        """Turn used for the rule and the cost: axis turn for a measured child, else the step's angle."""
        return _angle(din, self._out_dir(c, din, step)) if self.f["measured"][c] else _angle(din, step)

    def _passes(self, c, din, bound, step, first=False):
        """None if c passes thickness and direction from a fragment with direction din and bound, else the rule."""
        self.detail = None
        if not first and self.lower[c] > bound * (1 + TOL) + ABS_TOL:
            self.detail = (float(self.lower[c]), float(bound))
            return "thickness"
        if din is None:
            return None
        a = _angle(din, step)
        if a > STEP_MAX:
            self.detail = (a, np.nan)
            return "direction (backward)"
        if self.f["measured"][c]:
            t = self._turn(c, din, step)
            if t > TURN_MAX:
                self.detail = (t, np.nan)
                return "direction (bend)"
        return None

    def _coherent(self, s, c, din, bound, depth, seen):
        if depth == 0:
            return True
        f = self.f
        opts = []
        for g in self.nb[c]:
            if g in seen or self.is_seed[g] or (self.rct_owner[g] >= 0 and self.rct_owner[g] != s):
                continue
            o = self.owner[g]
            if o >= 0 and o != s:
                continue
            step = _unit(f["c"][g] - f["c"][c])
            if self._passes(g, din, bound, step) is None:
                opts.append((self._turn(g, din, step), g, step))
        for _, g, step in sorted(opts, key=lambda t: t[0]):       # look-ahead keeps the bound it started with
            if self._coherent(s, g, self._out_dir(g, din, step), bound, depth - 1, seen | {g}):
                return True
        return False

    def _new_state(self, hist, c, bound):
        h = hist + ((float(self.f["width"][c]),) if self.f["measured"][c] else ())
        h = h[-BOUND_K:]
        return h, (max(h) if h else bound)

    def _tie_check(self, s, c, cost):
        o = self.owner[c]
        if o >= 0 and o != s and cost < self.late.get((s, c), np.inf):
            self.late[(s, c)] = cost
        if o >= 0 and o != s and cost <= self.cost[c] + max(TIE_REL * self.cost[c], TIE_ABS):
            self.tie[c] = True
            self.tie_with[c] = s

    def run(self):
        f = self.f
        for s, fr in self.seeds.items():
            for c in fr:
                if self.owner[c] < 0:
                    self.owner[c], self.cost[c], self.is_seed[c] = s, 0.0, True
                    self.bound[c] = self.stem_diam[s]
        # stem -> first root fragments
        for s, fr in self.seeds.items():
            fr = [c for c in fr if self.owner[c] == s]
            for p in fr:
                for c in self.nb[p]:
                    if self.is_seed[c] or not self._ok_cand(s, c):
                        continue
                    step = _unit(f["c"][c] - f["c"][p])
                    rule = self._passes(c, None, self.stem_diam[s], step, first=True)
                    dout = self._out_dir(c, None, step)
                    if rule is None and not self._coherent(s, c, dout, self.stem_diam[s], LOOK, {p, c}):
                        rule, self.detail = "incoherent (first root)", None
                    if rule:
                        self.stop.setdefault((s, c), (rule, p) + (self.detail or (np.nan, np.nan)))
                        continue
                    hist, b = self._new_state((), c, self.stem_diam[s])
                    cost = float(np.linalg.norm(f["c"][c] - f["c"][p]))
                    heapq.heappush(self.heap, (cost, s, c, p, hist, b, tuple(dout)))
        while self.heap:
            cost, s, c, p, hist, b, din = heapq.heappop(self.heap)
            o = self.owner[c]
            if o >= 0:
                if o != s:
                    self._tie_check(s, c, cost)
                continue
            self.owner[c], self.cost[c], self.parent[c], self.bound[c] = s, cost, p, b
            self._expand(s, c, cost, hist, b, np.array(din))
        return self

    def _expand(self, s, p, cost, hist, bound, din):
        f = self.f
        cands = []
        for c in self.nb[p]:
            if self.owner[c] == s or self.is_seed[c] or not self._ok_cand(s, c):
                continue
            cands.append(c)
        if not cands:
            return
        steps = {c: _unit(f["c"][c] - f["c"][p]) for c in cands}
        short = [c for c in cands if not f["measured"][c]]
        fan = len(short) >= FAN_N and np.linalg.norm(np.mean([steps[c] for c in cands], 0)) < FAN_R
        ok = []
        for c in cands:
            rule = self._passes(c, din, bound, steps[c])
            detail = self.detail
            if (rule and rule.startswith("direction") and BRANCH_MAX is not None and not fan and f["measured"][c]
                    and max(_angle(din, steps[c]), self._turn(c, din, steps[c])) <= BRANCH_MAX
                    and self._coherent(s, c, self._out_dir(c, din, steps[c]), bound, LOOK, {p, c})):
                rule = None                              # branch
            if rule:
                self.stop.setdefault((s, c), (rule, p) + detail)
            else:
                ok.append((self._turn(c, din, steps[c]), c))
        ok.sort()
        acc = []
        if fan:
            for _, c in ok:
                if not acc and self._coherent(s, c, self._out_dir(c, din, steps[c]), bound, LOOK, {p, c}):
                    acc.append(c)
                else:
                    self.stop.setdefault((s, c), ("mat (fan)", p, len(cands), np.nan))
        else:
            tot = 0.0
            for _, c in ok:
                w2 = f["width"][c] ** 2 if f["measured"][c] else 0.0
                if acc and tot + w2 > (bound * (1 + SPLIT_TOL)) ** 2:
                    self.stop.setdefault((s, c), ("split sum", p, float(f["width"][c]), float(bound)))
                    continue
                tot += w2
                acc.append(c)
        for c in acc:
            turn = self._turn(c, din, steps[c])
            step_len = float(np.linalg.norm(f["c"][c] - f["c"][p]))
            cc = cost + step_len * (1 + (turn / TURN_MAX) ** 2)
            if self.owner[c] >= 0:                       # settled by another tree: only a tie check
                self._tie_check(s, c, cc)
                continue
            h, b = self._new_state(hist, c, bound)
            heapq.heappush(self.heap, (cc, s, c, p, h, b, tuple(self._out_dir(c, din, steps[c]))))

    def final_owner(self):
        """Owner per fragment; near-ties and everything grown from them -> -1 (flag in self.tie / tie_sub)."""
        own = self.owner.copy()
        order = np.argsort(self.cost)
        tie_sub = self.tie.copy()
        for c in order:
            p = self.parent[c]
            if p >= 0 and tie_sub[p] and self.owner[p] == self.owner[c]:
                tie_sub[c] = True
        own[tie_sub] = -1
        self.tie_sub = tie_sub
        return own


def thick_frags(gr):
    """Fragments where some tree's growth stopped on the thickness rule (candidates for touching roots)."""
    return np.array(sorted({c for (s, c), v in gr.stop.items() if v[0] == "thickness"}), np.int64)
