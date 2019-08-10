m = Healpix.Map{Int8, Healpix.RingOrder}(1)
for i in 1:length(m)
    m[i] = i
end

for i in m
    @test m[i] == i
end

@test m[1:3] == [1, 2, 3]
@test m[end-2:end] == [10, 11, 12]
