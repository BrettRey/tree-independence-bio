"""Independence polynomial computation and unimodality check for trees."""

try:
    from numpy import convolve as _np_convolve

    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

# int64 safe threshold: each output coeff is a sum of at most
# min(len(a),len(b)) products, each bounded by max_a * max_b.
# Keep total below 2^62 to avoid silent overflow in numpy int64.
_INT64_SAFE = 2**62


def _polymul_python(a: list[int], b: list[int]) -> list[int]:
    """Pure Python convolution with arbitrary-precision integers."""
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _polymul(a: list[int], b: list[int]) -> list[int]:
    if not a or not b:
        return []
    if _HAS_NUMPY:
        max_a = max(abs(x) for x in a)
        max_b = max(abs(x) for x in b)
        max_terms = min(len(a), len(b))
        if max_a > 0 and max_b > 0 and max_terms * max_a * max_b < _INT64_SAFE:
            return _np_convolve(a, b).astype(int).tolist()
    return _polymul_python(a, b)


def independence_poly(n: int, adj: list[list[int]]) -> list[int]:
    """Compute the independence polynomial of a tree.

    Uses iterative post-order DP to avoid recursion limits.

    Parameters
    ----------
    n : int
        Number of vertices.
    adj : list[list[int]]
        Adjacency list.

    Returns
    -------
    Coefficient list [i_0, i_1, ..., i_alpha] where i_k = number of
    independent sets of size k.
    """
    if n == 0:
        return [1]
    if n == 1:
        return [1, 1]

    # Handle forest by multiplying polynomials of components
    components_polys = []
    global_visited = [False] * n
    
    for start_node in range(n):
        if global_visited[start_node]:
            continue
            
        # Extract component
        component_nodes = []
        q = [start_node]
        global_visited[start_node] = True
        idx = 0
        while idx < len(q):
            u = q[idx]
            idx += 1
            component_nodes.append(u)
            for v in adj[u]:
                if not global_visited[v]:
                    global_visited[v] = True
                    q.append(v)
        
        # If single node component
        if len(component_nodes) == 1:
            components_polys.append([1, 1])
            continue
            
        # Relabel nodes for this component to 0..k-1
        k = len(component_nodes)
        mapping = {old: new for new, old in enumerate(component_nodes)}
        comp_adj = [[] for _ in range(k)]
        for u in component_nodes:
            for v in adj[u]:
                comp_adj[mapping[u]].append(mapping[v])
                
        # Compute poly for this component (it's a tree)
        # We can use the existing tree logic inline or recursive call
        # Let's use the tree logic here to avoid overhead of creating new list
        
        # Root at 0 (which is mapping[start_node])
        parent = [-1] * k
        children = [[] for _ in range(k)]
        visited = [False] * k
        order = []
        
        visited[0] = True
        bfs_queue = [0]
        head = 0
        while head < len(bfs_queue):
            v = bfs_queue[head]
            head += 1
            for u in comp_adj[v]:
                if not visited[u]:
                    visited[u] = True
                    parent[u] = v
                    children[v].append(u)
                    bfs_queue.append(u)
                    
        stack = [(0, False)]
        while stack:
            v, processed = stack.pop()
            if processed:
                order.append(v)
                continue
            stack.append((v, True))
            for c in children[v]:
                stack.append((c, False))
                
        dp0 = [[] for _ in range(k)]
        dp1 = [[] for _ in range(k)]
        
        for v in order:
            if not children[v]:
                dp0[v] = [1]
                dp1[v] = [0, 1]
            else:
                prod = [1]
                for c in children[v]:
                    summand = _polyadd(dp0[c], dp1[c])
                    prod = _polymul(prod, summand)
                dp0[v] = prod
                
                prod = [1]
                for c in children[v]:
                    prod = _polymul(prod, dp0[c])
                dp1[v] = [0] + prod
                
        components_polys.append(_polyadd(dp0[0], dp1[0]))

    # Multiply all component polynomials
    if not components_polys:
        return [1]
    
    result = components_polys[0]
    for p in components_polys[1:]:
        result = _polymul(result, p)
        
    return result


def _polyadd(a: list[int], b: list[int]) -> list[int]:
    """Add two polynomials represented as coefficient lists."""
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def is_unimodal(seq: list[int]) -> bool:
    """Check whether a sequence is unimodal.

    A sequence is unimodal if it never increases after a decrease:
    no index i with seq[i-1] > seq[i] < seq[i+1].

    More precisely: once the sequence starts decreasing, it never
    increases again.
    """
    if len(seq) <= 2:
        return True
    decreasing = False
    for i in range(1, len(seq)):
        if seq[i] > seq[i - 1]:
            if decreasing:
                return False
        elif seq[i] < seq[i - 1]:
            decreasing = True
    return True


def is_log_concave(seq: list[int]) -> bool:
    """Check whether a sequence is log-concave.

    A sequence is log-concave if a_k^2 >= a_{k-1} * a_{k+1} for all
    valid k. Uses integer arithmetic to avoid floating-point issues.
    """
    for k in range(1, len(seq) - 1):
        if seq[k] * seq[k] < seq[k - 1] * seq[k + 1]:
            return False
    return True


def log_concavity_ratio(seq: list[int]) -> tuple[float, int]:
    """Measure worst log-concavity violation.

    Returns (max_ratio, position) where ratio = a_{k-1}*a_{k+1} / a_k^2.
    Ratio > 1 means log-concavity fails at that position.
    Larger ratio = worse violation.
    """
    worst = 0.0
    worst_k = -1
    for k in range(1, len(seq) - 1):
        sq = seq[k] * seq[k]
        if sq == 0:
            continue
        ratio = (seq[k - 1] * seq[k + 1]) / sq
        if ratio > worst:
            worst = ratio
            worst_k = k
    return worst, worst_k


def near_miss_ratio(seq: list[int]) -> tuple[float, int]:
    """Measure how close a unimodal sequence is to violating unimodality.

    After the first strict decrease, tracks the maximum ratio
    a_{j+1}/a_j in the "should be non-increasing" tail. A ratio > 1
    would be a unimodality violation. The closer to 1, the closer to
    a counterexample.

    Returns (max_ratio, position_j) or (0.0, -1) if the sequence
    never decreases.
    """
    if len(seq) <= 2:
        return 0.0, -1

    # Find first strict decrease
    first_descent = -1
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            first_descent = i
            break

    if first_descent == -1:
        return 0.0, -1  # monotone non-decreasing

    # In the descending tail, find max a_{j+1}/a_j
    worst = 0.0
    worst_j = -1
    for j in range(first_descent, len(seq) - 1):
        if seq[j] == 0:
            continue
        ratio = seq[j + 1] / seq[j]
        if ratio > worst:
            worst = ratio
            worst_j = j
    return worst, worst_j
