# TE11CylCavity - Find ϵᵣ and Q for dielectric sample inserted in a TE₁₁ₚ cavity

## Summary as of 5/30/2025: Comparing to CST for Case2 Geometry
For case2, the best results for obtaining the actual CST dielectric constant and loss tangent are
obtained using test_findcase2fresQ.jl, compared to test2_findcase2fresQ.jl and test3_findcase2fresQ.jl.
For the best method the conductivity and cavity length are first adjusted slightly so that the S-parameter
model exactly agrees with the CST prediction for the empty cavity.  These new values are then used when searching
for the dielectric properties.

test3_findcase2fresQ.jl uses a similar technique, but it adjust a (cylinder radius) rather than d (length). Its results
are poorer for dielectric constant.

test2_findcase2fresQ.jl doesn't adjust the dimensions or conductivity, but it does adjust the frequency and Q goals.  This
didn't work very well compared to the other two techniques.

# Summary as of 6/1/2025
- Added interpolation of rectangular/cylindrical waveguide transition S-parameters.
- Retained only scripts test_findcase2afresQ.jl and test_findcase3fresQ.jl.  These use more accurate CST analyses
  from Mike.  They also use only the technique of adjusting conductivity and cavity length `d` to make the S-parameter
  analysis of the empty cavity agree with the given measured or computed resonant frequency and Q.
- The predicted values for Case2a are `(ϵᵣ = 5.75014044812833, tanδ = 0.0032051587128061635` while the 
  actual values are (ϵᵣ = 5.7, tanδ = 0.003).
- Mike sent me Case3 without informing me what the actual values are.  Running the script produces the prediction
  `(ϵᵣ = 7.979918953136148, tanδ = 0.003990828375213885`. After I reported these to Mike he informed me that the 
  actuals are `(ϵᵣ = 8.0, tanδ = 0.004)`.  Not bad!  The better accuracy of the Case3 prediction is likely due to the
  fact that Mike centered the sample right at one of the peaks of the TE117 cavity modes.  For Case2a, the sample is 
  located pretty near the end of the cavity where the transverse E-field vanishes (approximately, due to the coupling iris
  and the finite conductivity of the shorting wall).

