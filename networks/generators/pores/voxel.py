import numpy as np


class VoxelGrid:
    def __init__(self, bounds:tuple[float, float, float], dmin:float):
        self.bounds = bounds

        self.rmin = dmin/2
        max_voxel_size = np.sqrt(dmin**2/3) # *(1-1e-4)
        
        self.lengths = (
            int(np.ceil(bounds[0] / max_voxel_size)),
            int(np.ceil(bounds[1] / max_voxel_size)),
            int(np.ceil(bounds[2] / max_voxel_size)),
        )
        self.voxel_sizes = (
            bounds[0] / self.lengths[0],
            bounds[1] / self.lengths[1],
            bounds[2] / self.lengths[2]
        )

        self.min_voxel_size = min(self.voxel_sizes)

        self.length = self.lengths[0] * self.lengths[1] * self.lengths[2]
        self.indices = set()
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
            (index[0]+0.5) * self.voxel_sizes[0], 
            (index[1]+0.5) * self.voxel_sizes[1], 
            (index[2]+0.5) * self.voxel_sizes[2]
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
            index[0] * self.voxel_sizes[0], 
            index[1] * self.voxel_sizes[1], 
            index[2] * self.voxel_sizes[2]
        )
        bounds_max = (
            (index[0]+1) * self.voxel_sizes[0], 
            (index[1]+1) * self.voxel_sizes[1], 
            (index[2]+1) * self.voxel_sizes[2]
        )
        return (bounds_min, bounds_max)

    def get_index(self, pos: tuple[float, float, float]) -> tuple[float, float, float]:
        index = (
            int(np.floor(pos[0]/self.voxel_sizes[0])),
            int(np.floor(pos[1]/self.voxel_sizes[1])),
            int(np.floor(pos[2]/self.voxel_sizes[2])),
        )
        return index


    def remove_filled_voxels(
        self,
        center: tuple[float, float, float],
        radius: float,
        ) -> set[tuple[int, int, int]]:
        voxel_radius = int(np.ceil(radius/self.min_voxel_size)) + 1
        index_center = self.get_index(center)
        radius_squared = radius**2
        assert index_center in self.indices
        removed_indices = set()
        # print("remove ", len_start := len(self.indices), " | ", end="")

        for i in range(2* voxel_radius):
            index_x = index_center[0] + i - voxel_radius
            if index_x < 0 or index_x >= self.lengths[0]:
                continue
            for j in range(2* voxel_radius):
                index_y = index_center[1] + j - voxel_radius
                if index_y < 0 or index_y >= self.lengths[1]:
                    continue
                for k in range(2* voxel_radius):
                    index_z = index_center[2] + k - voxel_radius
                    if index_z < 0 or index_z >= self.lengths[2]:
                        continue
                    index_voxel = (
                        index_x,
                        index_y,
                        index_z
                    )
                    center_voxel = self.get_voxel_center(index_voxel)
                    distance_squared = (
                        (center_voxel[0] - center[0])**2 +
                        (center_voxel[1] - center[1])**2 +
                        (center_voxel[2] - center[2])**2
                    ) 
                    if distance_squared < radius_squared:
                        removed_indices.add(index_voxel)
                        self.indices.discard(index_voxel)

        # print(len_end := len(self.indices), f" (deleted: {len_start-len_end})")
        return removed_indices