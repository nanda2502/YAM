{-# LANGUAGE FlexibleContexts #-}
module Utils (
    writeMatrixToCSV,
    formatResult
) where

import Data.List (intercalate)
import Text.Printf (printf)

-- | Writes the matrix to a CSV file
writeMatrixToCSV :: FilePath -> [[Double]] -> IO ()
writeMatrixToCSV fileName matrix = do
    let csvData = map (intercalate "," . map (printf "%.2f")) matrix
    writeFile fileName $ unlines csvData

-- | Helper function to format the result into a CSV line with rounding
formatResult :: (Int, String, Double, String, Double) -> String
formatResult (n, adjMatrixBinary, alpha, strategy, expectedSteps) =
    intercalate "," [show n, adjMatrixBinary, show alpha, strategy, printf "%.2f" expectedSteps]