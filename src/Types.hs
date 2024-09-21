{-# LANGUAGE FlexibleContexts #-}
module Types (
    Trait,
    Repertoire,
    PayoffVector,
    AdjacencyMatrix,
    Strategy(..)
) where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V

type Trait = Int
type Repertoire = VU.Vector Bool
type PayoffVector = VU.Vector Double
type AdjacencyMatrix = V.Vector (VU.Vector Bool)

-- | Learning strategies
data Strategy = RandomLearning | PayoffBasedLearning deriving (Eq, Show)