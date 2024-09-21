{-# LANGUAGE FlexibleContexts #-}
module LinearAlgebra (
    solveLinearSystem
) where

import qualified Data.Matrix as Mat
import qualified Data.Vector as V

-- | Solves a linear system of equations using matrix inversion
solveLinearSystem :: [[Double]] -> [Double] -> Maybe [Double]
solveLinearSystem a b =
    let aMatrix = Mat.fromLists a                     -- Convert the list of lists to a Matrix
        bVector = V.fromList b                        -- Convert the list to a Data.Vector
        bMatrix = Mat.colVector bVector               -- Convert the Vector to a column Matrix
    in case Mat.inverse aMatrix of
        Left _ -> Nothing                             -- The matrix is singular or not invertible
        Right invA ->
            let xMatrix = Mat.multStd invA bMatrix    -- Matrix multiplication using multStd
                xList = concat $ Mat.toLists xMatrix  -- Convert the solution back to a list
            in Just xList