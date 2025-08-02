class KDNode:
    def __init__(self, point, left=None, right=None):
        self.point = point
        self.left = left
        self.right = right

def build_kdtree(points, depth=0):
    n = len(points)
    if n == 0:
        return None
    k = points.shape[1]  # Assuming all points have the same dimension
    axis = depth % k

    # Median of points sorted by the current axis
    indices = np.argsort(points[:, axis])
    median_idx = n // 2

    # Create node and construct subtrees using indices to avoid data copy
    return KDNode(
        point=points[indices[median_idx]],
        left=build_kdtree(points[indices[:median_idx]], depth + 1),
        right=build_kdtree(points[indices[median_idx + 1:]], depth + 1)
    )

def distance_squared(point1, point2):
    return np.sum((np.array(point1) - np.array(point2)) ** 2)

def query_kdtree(root, point, radius, depth=0):
    if root is None:
        return []
    k = len(point)
    axis = depth % k
    points = []
    
    # Check the distance to the point in this node
    if np.sum((np.array(root.point) - np.array(point)) ** 2) <= radius**2:
        points.append(root.point)
    
    # Determine which branch to search
    next_branch = None
    opposite_branch = None
    if point[axis] < root.point[axis]:
        next_branch = root.left
        opposite_branch = root.right
    else:
        next_branch = root.right
        opposite_branch = root.left
    
    # Search the next branch
    points.extend(query_kdtree(next_branch, point, radius, depth + 1))
    
    # Check whether we should search the opposite branch
    if (point[axis] - root.point[axis])**2 <= radius**2:
        points.extend(query_kdtree(opposite_branch, point, radius, depth + 1))
    
    return points