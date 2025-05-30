# Summary as of 5/30/2025
For case2, the best results for obtaining the actual CST dielectric constant and loss tangent are
obtained using test_findcase2fresQ.jl, compared to test2_findcase2fresQ.jl and test3_findcase2fresQ.jl.
For the best method the conductivity and cavity length are first adjusted slightly so that the S-parameter
model exactly agrees with the CST prediction for the empty cavity.  These new values are then used when searching
for the dielectric properties.

test3_findcase2fresQ.jl uses a similar technique, but it adjust a (cylinder radius) rather than d (length). Its results
are poorer for dielectric constant.

test2_findcase2fresQ.jl doesn't adjust the dimensions or conductivity, but it does adjust the frequency and Q goals.  This
didn't work very well compared to the other two techniques.


