import unittest
import random
import simpleMap as sM
import numpy

"""The following unitests check the correctness of implementations of 
the functions to be completed.
"""

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.testNo = 100
    
    def testMinimizerIndex(self):
        """Tests cA.Minimizer indexer using a couple of hand crafted examples.
        """
        w, k = 5, 3
        
        def getMinmers(mI):
            minmersFound = set()
            for minmer in mI.minimizerMap:
                for minmerIndex in mI.minimizerMap[minmer]:
                    minmersFound.add((minmer, minmerIndex))
            return minmersFound
        
        x = "TACCCCTCAGATGCTTAAGC"
        xM = set([('ACC', 1), ('AAG', 16), ('CCC', 3), ('CCC', 2), ('CTT', 13), ('CCT', 4), ('ATG', 10), ('CAG', 7), ('AGA', 8)])
        mI = sM.MinimizerIndexer(x, w=w, k=k, t=1000) # Set t so that no filtering takes place
        # Check the sets of minimizer sets are equal        
        self.assertEqual(xM, getMinmers(mI))
        
        x = "TTCGTTTTTCACTCCGGGTC"
        xM = set([('CAC', 9), ('TCA', 8), ('GTT', 3), ('TTC', 7), ('GGG', 15), ('ACT', 10), ('CGT', 2), ('CGG', 14), ('TTT', 4), ('CCG', 13)])
        mI = sM.MinimizerIndexer(x, w=w, k=k, t=1000) # Set t so that no filtering takes place
        # Check the sets of minimizer sets are equal        
        self.assertEqual(xM, getMinmers(mI))

    def testSeedCoordinates(self):
        """Tests cA.SeedCluster indexer using randomly generated 
        examples, checking the output using a simple brute force algorithm.
        """ 
        # Ensure this doesn't pass easily without being correct
        for test in range(self.testNo): #Do self.testNo random tests
            xLen = random.choice(range(0, 10))
            yLen = random.choice(range(0, 10))
            l = random.choice(range(0, 10))
            
            # Make random set of seeds
            seedList = []
            seedComponents = {}
            for x in range(xLen):
                if random.random() > 0.5:
                    ys = []
                    for y in range(yLen):
                        if random.random() > 0.9:
                            ys.append(y)
                            seedComponents[(x, y)] = [(x, y)]
                    seedList.append((x, tuple(ys)))
            
            #print "l", l, "xLen", xLen, "yLen", yLen, "seedList", seedList
            
            # Run stupid brute force algorithm to collate seeds
            # Note this algorithm is very expensive scaling roughly with the fourth power 
            # of the length of the input sequences
            seeds = list(seedComponents.keys())
            for i in range(len(seeds)):
                x1, y1 = seeds[i]
                for j in range(i+1, len(seeds)):
                    x2, y2 = seeds[j]
                    if abs(x1 - x2) <= l and abs(y1 - y2) <= l:
                        component1 = seedComponents[(x1, y1)]
                        component2 = seedComponents[(x2, y2)]
                        # Merge components if not already merged
                        if component1 != component2:
                            component1 += component2
                            for (x3, y3) in component2:
                                seedComponents[(x3, y3)] = component1
            
            # Generate seed clusters
            seedClusters = sM.SeedCluster.clusterSeeds(seedList, l)
            
            # Check we get what we expect
            for seedCluster in seedClusters:
                expectedCluster = list(seedComponents[seedCluster.seeds[0]])
                expectedCluster.sort()
                self.assertEquals(expectedCluster, seedCluster.seeds)            
      
    def testSmithWaterman(self):
        """Tests cA.SmithWaterman using randomly generated 
        examples, checking the output using a simple brute force algorithm.
        """ 
        gapScore=-2
        matchScore=3
        mismatchScore=-3
        
        x = "ATTCTGC"
        y = "CGTTGCCATCTCCGTATGA"
        maxScore = 11
        optimalAlignment = [(0, 7), (2, 8), (3, 9), (4, 10), (6, 11)]
        
        # Get alignment object
        alignment = sM.SmithWaterman(x, y, gapScore=gapScore, matchScore=matchScore, mismatchScore=mismatchScore)
        # Check that alignment scores agree
        self.assertEqual(maxScore, alignment.getMaxAlignmentScore())
        # Check that tracebacks agree
        self.assertEqual(optimalAlignment, alignment.getAlignment())
        
        x = "AGAAGCGAGATTACAC"
        y = "CCT"
        maxScore = 4
        optimalAlignment = [(13, 0), (15, 1)]
        
        # Get alignment object
        alignment = sM.SmithWaterman(x, y, gapScore=gapScore, matchScore=matchScore, mismatchScore=mismatchScore)
        # Check that alignment scores agree
        self.assertEqual(maxScore, alignment.getMaxAlignmentScore())
        # Check that tracebacks agree
        self.assertEqual(optimalAlignment, alignment.getAlignment())


if __name__ == '__main__':
    unittest.main()