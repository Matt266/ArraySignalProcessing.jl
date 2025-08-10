struct AzEl <: PlaneWave
    coords::AbstractMatrix
    function AzEl(coords)
         # coords = [azimtuh, elevation] 2xD matrix

        if ndims(coords) == 0
            coords = [coords; 0]
        end

        # only when single az/el pair
        # transpose changes ndims 
        # az-only: [az0, az2]', az/el-pair: [az, el]
        # ndims(transpose([3, 4])) is 2
        # ndims([3, 4]) is 1
        if ndims(coords) == 1
            coords = reshape(coords, 2, 1)
        end

        M, D = size(coords)

        if M == 1
            coords = [coords; zeros(1, D)]
        end

        return new(coords)
    end
end