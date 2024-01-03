import numpy as np


class VoxelGrid:
    def __init__(self, bounds:tuple[float, float, float], dmin:float):
        self.bounds = bounds

        self.voxel_size = np.sqrt(dmin**2/3)
        self.lengths = (
            int(np.ceil(bounds[0] / self.voxel_size)),
            int(np.ceil(bounds[1] / self.voxel_size)),
            int(np.ceil(bounds[2] / self.voxel_size)),
        )

        self.length = self.lengths[0] * self.lengths[1] * self.lengths[2]
        self.indices = {}
        for i in range(self.lengths[0]):
            for j in range(self.lengths[1]):
                for k in range(self.lengths[2]):
                    self.indices.add((i, j, k))

    def get_voxel_center(
        self, 
        index: tuple[int, int, int]
        )->tuple[float, float, float]:
        for i in range(3):
            assert index[i] >= 0 and index[i] < self.lengths[i]
        center = (
            (index[0]+0.5) * self.voxel_size, 
            (index[1]+0.5) * self.voxel_size, 
            (index[2]+0.5) * self.voxel_size
        )
        return center
        

    def get_voxel_bounds(
        self, 
        index: tuple[int, int, int]
        )->tuple[
            tuple[float, float, float], 
            tuple[float, float, float]
        ]:
        for i in range(3):
            assert index[i] >= 0 and index[i] < self.lengths[i]

        bounds_min = (
            index[0] * self.voxel_size, 
            index[1] * self.voxel_size, 
            index[2] * self.voxel_size
        )
        bounds_max = (
            (index[0]+1) * self.voxel_size, 
            (index[1]+1) * self.voxel_size, 
            (index[2]+1) * self.voxel_size
        )
        return (bounds_min, bounds_max)

    def get_index(self, pos: tuple[float, float, float]) -> tuple[float, float, float]:
        index = (
            int(np.floor(pos[0]/self.voxel_size)),
            int(np.floor(pos[1]/self.voxel_size)),
            int(np.floor(pos[2]/self.voxel_size)),
        )
        return index


    def remove_filled_voxels(
        self,
        center: tuple[float, float, float],
        radius: float,
        ):
        voxel_radius = int(np.ceil(radius/self.voxel_size))
        index_center = self.get_index(center)
        radius_squared = radius**2

        for i in range(2* voxel_radius + 1):
            if index_center[0] + i - voxel_radius < 0:
                continue
            for j in range(2* voxel_radius + 1): 
                if index_center[1] + j - voxel_radius < 0:
                    continue
                for k in range(2* voxel_radius + 1):
                    if index_center[2] + k - voxel_radius < 0:
                        continue
                    index_voxel = (
                        index_center[0] + i - voxel_radius,
                        index_center[1] + j - voxel_radius,
                        index_center[2] + k - voxel_radius
                    )
                    center_voxel = self.get_voxel_center(index_voxel)
                    distance_squared = (
                        (center_voxel[0] - center[0])**2 +
                        (center_voxel[1] - center[0])**2 +
                        (center_voxel[2] - center[0])**2
                    )
                    if distance_squared < radius_squared:
                        self.indices.remove(index_voxel)