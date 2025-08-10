struct WaveVec <: PlaneWave
    coords::AbstractMatrix
    function WaveVec(coords)
        # coords = [kx; ky; kz] 3xD matrix

        if ndims(coords) == 0
            coords = [coords; 0; 0]
        end

        # only when single k vector
        # transpose changes ndims
        # single k vector: [kx, ky, kz]'
        # multiple kx: [kx, kx, kx]
        if ndims(coords) == 1
            coords = reshape(coords, 3, 1)
        end

        M, D = size(coords)

        if M == 1
            coords = [coords; zeros(2, D)]
        elseif M==2
            coords = [coords; zeros(1, D)]
        end 

        return new(coords)
    end
end