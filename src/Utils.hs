{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE LambdaCase #-}
module Utils (
    writeMatrixToCSV,
    formatResult,
    readAdjacencyMatrices
) where

import Data.List (intercalate)
import Text.Printf (printf)
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import Types (AdjacencyMatrix)
import Control.Monad (forM_)
import Data.Either (rights)

-- | Writes the matrix to a CSV file
writeMatrixToCSV :: FilePath -> [[Double]] -> IO ()
writeMatrixToCSV fileName matrix = do
    let csvData = map (intercalate "," . map (printf "%.2f")) matrix
    writeFile fileName $ unlines csvData

-- | Helper function to format the result into a CSV line with rounding
formatResult :: (Int, String, Double, String, Double) -> String
formatResult (n, adjMatrixBinary, alpha, strategy, expectedSteps) =
    intercalate "," [show n, adjMatrixBinary, show alpha, strategy, printf "%.2f" expectedSteps]

-- | Reads adjacency matrices from a CSV file where each line is a binary string
readAdjacencyMatrices :: Int -> IO [AdjacencyMatrix]
readAdjacencyMatrices n = do
    let filePath = "./data/adj_mat_" ++ show n ++ ".csv"
    content <- readFile filePath
    let linesList = lines content
    let matrices = map (binaryStringToAdjMatrix n) linesList
    forM_ matrices (\case
        Left err -> error err
        Right _ -> return ())
    return (rights matrices)

-- | Converts a binary string to an adjacency matrix
binaryStringToAdjMatrix :: Int -> String -> Either String AdjacencyMatrix
binaryStringToAdjMatrix n str =
    let bits = mapM charToBool (filter (/= ',') str)
        totalBits = n * n
    in case bits of
        Left err -> Left err
        Right bitList ->
            if length bitList /= totalBits
                then Left $ "Invalid length: Expected " ++ show totalBits ++ " bits, got " ++ show (length bitList)
                else Right $ V.fromList [ VU.fromList (take n (drop (i * n) bitList)) | i <- [0 .. n - 1] ]

-- | Helper function to convert a character to a Bool
charToBool :: Char -> Either String Bool
charToBool '1' = Right True
charToBool '0' = Right False
charToBool c   = Left $ "Invalid character in binary string: " ++ [c]