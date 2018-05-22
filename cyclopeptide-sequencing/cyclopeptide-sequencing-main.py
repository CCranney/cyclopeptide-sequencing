import sys, operator

infile = str(sys.argv[1])
charge = int(sys.argv[2])

with open( infile ) as file:
#with open( 'unknownExperimentalSpectrum.txt' ) as file:
    data = [ float(i.strip())-charge for i in file ]


###################################
#
#
#   Section 1:
#   Global Variables 
#
#
##################################

# AminoAcids_String_Weight is a dictionary with the following contents:
#
# Keys are letters representing amino acids
# Values are corresponding weights in daltons
AminoAcids_String_Weight = dict()
with open( "aminoAcidWeights.txt" ) as file:
    AminoAcids_String_Weight = dict( [line.split() for line in file] )
# Converts values from strings to floats
for i in AminoAcids_String_Weight:
    AminoAcids_String_Weight[i] = float(AminoAcids_String_Weight[i])

##################################
#
#
#   Section 2:
#   Basic Functions Used
#   Throughout the Program
#
#
##################################

# Determines if two float values are within "tolerance" of each other.
#
# x is first float
# y is second float
# Natural default tolerance is 0.3
#   Returns True if approximately equal,
#   Returns False otherwise

def Approx( x, y, tolerance=0.300000001 ):
    return abs(x-y) <= tolerance

# Runs Approx on a float and each element in a set.
#
# value is the float
# setCompare is the set value is being compared to
# Natural default tolerance is 0.3
#   If the value in a setCompare matches, returns matching set value
#   If not, returns -1
def ApproxSet( value, setCompare, tolerance=0.300000001 ):
    for i in setCompare:
        if Approx( value, i ,tolerance = tolerance):
            return i
    return -1

# Given a list of floats or integers, this function will group
# similarly sized items together and return their average,
# thereby eliminating approximate duplicates from the list.
#
# listWithDuplicates is list of floats or ints that may have duplicates
#   Returns a list of tuples, with the first item being the 
#       number of times the approximate value appeared in the 
#       original listWithDuplicates.
# 
#Note: tolerance for set ApproxSet() set to 0.6. The variation of each
# data point is not accepted beyond 0.3, so by adding two together
# their combined variation would not exceed 0.6.
def EliminateApproximateDuplicateValues( listWithDuplicates ):
    d = dict()
    for i in listWithDuplicates:
        m = ApproxSet( i, set(d.keys()), tolerance = 0.60000000001 )
        if m == -1:
            d[i] = 1
        else:
            tempM = m
            tempN = d[m]
            tempM *= tempN
            tempN += 1
            tempM += i
            tempM /= float(tempN)
            d.pop(m)
            d[tempM] = tempN
    d_score = sorted(d.items(), key=operator.itemgetter(1),reverse=True)
    return d_score

# This function takes all data points and returns the "opposite" segment.
#  If the data point is a true representative of the cyclopeptide, then full-dp
#  would also be a data point. If the data point is false, another false data
#  point would be created instead.
#
# refinedData is a set containing the experimental spectrum. So named because
#  some actions can be taken to eliminate false data points, such as removing
#  all cyclopeptides within 50 of the full sequence.
# full is the predicted weight of the full cyclopeptide.
# 
# returns a set containing all elements in refined Data and their counterparts
def DoubleTheData( refinedData, full ):
    doubleData = set()
    for i in refinedData:
        temp = ApproxSet( full-i, refinedData )
        doubleData.add(i)
        if temp == -1:
            doubleData.add( full-i )
    return doubleData

#NOTE I'm trying to get rid of this function, so I won't comment on it as of yet
def CreateSmallSegments( esData, full ):
    smallSegments = set()
    for i in esData:
        temp = full-i
        m = ApproxSet( i, smallSegments )
        n = ApproxSet( temp, smallSegments )
        if abs(temp - i) < 100 and m == -1 and n == -1:
            smallSegments.add(i)
            smallSegments.add(temp)
        elif temp > i and m == -1:
            smallSegments.add(i)
        elif temp < i and n == -1:
            smallSegments.add(temp)
    return smallSegments

# Returns the theoretical spectrum of a cyclopeptide chain. This is largely used to
#  test the validity of an incomplete cyclopeptide chain. Presumably, the closer the chain is
#  to representing the original cyclopeptide, the more it's own theoretical spectrum
#  will appear in the data set.
#  
# chain is the cyclopeptide chain being tested.
#
# returns ts, the theoretical spectrum of the chain.
def TheoreticalSpectrumChain( chain ):
    ts = []
    for i in xrange( len(chain) ):
        newPeptide = list(chain[i:]) + list(chain[0:i])
        for j in xrange( 1,len(newPeptide) ):
            ts.append(sum(newPeptide[0:j]))
    return ts

# Same as TheoreticalSpectrumChain, except calculated the theoretical spectrum of a linear
#  peptide, not a cyclopeptide
def LinearTheoreticalSpectrum( chain ):
    tsc = []
    for k in xrange(len(chain)):
        tsc += [ sum( chain[k:i] ) for i in xrange(k+1, len(chain)+1 ) ]
    return tsc

# Scores a chain by comparing it's theoretical spectrum to the experimental spectrum (data set).
#
# chain is the cyclopeptide chain being tested.
# esDataSet is the experimental spectrum.
#
# returns count, an int representing the number of times a segment of the theoretical spectrum
#  of the chain matched a value in the experimental spectrum.
#NOTE I'm considering setting the tolerance to 0.6. Do it further down the road of corrections, though.
def Score( chain, esDataSet ):
    esDataSetCopy = set( sorted( esDataSet ) )
    tsc = TheoreticalSpectrumChain( chain ) 
    count = 0
    for i in tsc:
        match = ApproxSet( i, esDataSetCopy )
        if match!= -1:
            esDataSetCopy.remove(match)
            count += 1
    return count


# Same as Score(), except this algorithm tests the chain as a LINEAR peptide, not a cyclopeptide.
#  The difference is important.
#
# chain is the linear peptide chain being tested.
# esDataSet is the experimental spectrum.
#
# returns count, just as in Score()
def ScoreLinear( chain, esDataSet ):
    esDataSetCopy = set( sorted( esDataSet ) )
    tsc = [ sum( chain[0:i] ) for i in xrange(1, len(chain) ) ]
    tsc += [ sum(chain[i:len(chain)]  ) for i in xrange( len(chain) )]
    count = 0
    for i in tsc:
        match = ApproxSet( i, esDataSetCopy )
        if match != -1:
            esDataSetCopy.remove(match)
            count += 1
    return count

# Takes a cyclopeptide chain and returns all of it's different possible orientations.
#  This includes the original, the same chain starting from each other segment, and the reverse
#  of each.
#
# chain represents a cyclopeptide sequence.
#
# returns matches, a list of all different possible orientations of chain.
def FindMatchingSegments( chain ):
    matches = [chain]
    for i in xrange( len(chain) ): 
        newChain = list(chain[i:]) + list(chain[0:i])
        reverseChain = newChain[::-1]
        matches.append(newChain)
        matches.append(reverseChain)
    return matches

# Tests whether two chains are approximately equal. This function assumes that the two chains are
#  already lined up.
#
# chain1 represents one chain of segments.
# chain2 represents a second chain of segments.
#
# Returns True if they match, False if otherwise
def AreChainsSame( chain1, chain2 ):
    if len(chain1) != len(chain2):
        return False
    for i in xrange( len(chain1) ): 
        if Approx( chain1[i],chain2[i], tolerance = 0.60000000000001 ) == False:
            return False
    return True

# Takes a list of chains and removes all duplicates from the list. This was slightly 
#  tougher than it sounds, as the same chain could appear in a different orientation.
#  The following is one methodology that works, but is incredibly slow if the list of
#  chains is too long. A one-on-one comparison method may work better, but this will do
#  in a pinch.
#
# chains is a list of chains that may have duplicates.
#
# returns c, a list of chains with the duplicates removed.
def EliminateDuplicateChains( chains ):
    c = []
    allChains = []
    for i in chains:
        check = True    
        for j in allChains:
            if AreChainsSame( i, j ) == True:
                check = False
        if check:
            allChains += FindMatchingSegments( i )
            c.append(i)
    return c

# This algorithm is largely taken from the Bioinformatics Algorithm book. It takes a list 
#  of viable cyclopeptides, orders them by their Score() function value, and then removes
#  all cyclopeptides who's score is below the cyclopeptide in the N'th place.
# 
# leaderboard is the list of cyclopeptide chains being compared.
# esData is the experimental spectrum data set.
# N is the cutoff number described in the comment above.
#
# returns a trimmed leaderboard.
def Trim( leaderboard, esData, N ):
    if N > len(leaderboard):
        return leaderboard
    scores = [ 0 for i in leaderboard ]
    for i in xrange( len(leaderboard) ):
        scores[i] = Score( leaderboard[i], set(esData) )
    scores, leaderboard = (list(x) for x in zip(*sorted(zip(scores, leaderboard), key=lambda pair: pair[0])))
    for i in xrange( len(scores)-1,-1,-1 ):
        if scores[i] < scores[len(scores)-N]:
            leaderboard = leaderboard[i+1:]
            break
    return leaderboard
   
# Same as Trim(), only this one utilizes the ScoreLinear() function instead of the Score() function.
#  In essence, this is used to Trim a list of LINEAR peptides, not cyclopeptides.
def TrimLinear( leaderboard, esData, N ):
    if N > len(leaderboard):
        return leaderboard
    scores = [ 0 for i in leaderboard ]
    for i in xrange( len(leaderboard) ):
        scores[i] = ScoreLinear( leaderboard[i], set(esData) )
    scores, leaderboard = (list(x) for x in zip(*sorted(zip(scores, leaderboard), key=lambda pair: pair[0])))
    for i in xrange( len(scores)-1,-1,-1 ):
        if scores[i] < scores[len(scores)-N]:
            leaderboard = leaderboard[i+1:]
            break
    return leaderboard

##################################
#
#
#   Section 3:
#   Predicting the Full Size
#    Of a Cyclopeptide given
#    an Experimental Spectrum
#    of Mass Spectrometry Data
#
#
##################################

# Knowing the full size of the cyclopeptide is crucial to the 
#  sequencing process. The following function predicts several
#  candidate full sequences and matching candidate amino acids
#  that are more likely to appear in the final sequence.
#
# esData is the experimental data set. This is in list format.
# N is a trimming cutoff number. In the first part, the top
#  N candidates will be kept, as well as any candidates who has
#  a score equal to the N'th candidate.
# variationFromMaximum is another trimming cutoff number. While
#  technically used before N, this one is a variable I am considering
#  doing away with. In essence, all possible candidates outside 
#  variationFromMaximum from the Max value in esData are eliminated
#  from consideration.
# 
# Returns a list of tuples.
#  The first item in the tuple is the candidate full length.
#  The second item in the tuple is a string of amino acids
#   that could correspond to the full length.
def PredictFullSizes( esData, N, variationFromMaximum ):
    # Part One: Generating list of candidate complete sequences.
    #  This is run on the calculation that adding two opposite fragments
    #  of the full cyclopeptide would result in the full peptide. 
    #  As such, adding each element of the experimental data set 
    #  would result in the full peptide weight appearing several times
    #  more than simply by random chance.
    
    # adding all fragments together and storing results in a list
    allPossibleFullSizes = []
    for i in xrange( len( esData ) ):
        for j in xrange( i+1, len(esData) ):
            allPossibleFullSizes.append(esData[i] + esData[j])
    
    # Calculating frequently appearing larger sequences
    sorted_score = EliminateApproximateDuplicateValues( allPossibleFullSizes )

    # Eliminating all candidate sequences outside of the range indicated by variationFromMaximum
    for i in xrange( len(sorted_score)-1,-1,-1 ):
        if sorted_score[i][0] < max( esData )-variationFromMaximum or sorted_score[i][0] > max(esData) + variationFromMaximum:
            sorted_score.pop(i)
    
    # Eliminating all candidate sequences with a score less than the N'th candidate
    for i in xrange( len(sorted_score) ):
        if sorted_score[i][1] < sorted_score[N][1]:
            sorted_score = sorted_score[0:i]
            break

    # Part Two: Generating candidate amino acids attached to the candidate sequences.
    #  This step serves three purposes: normalizing the candidate sequence, further
    #  eliminating false candidates, and providing likely amino acids that could be found
    #  in the cyclopeptide.

    # Attach matching amino acids to each candidate sequence.
    topScores = []
    for seq in sorted_score:
        temp = [ seq[0], '' ]
        for AA in AminoAcids_String_Weight:
            if ApproxSet( seq[0]-AminoAcids_String_Weight[AA], set(esData), tolerance = 1.000000001 ) != -1:
                temp[1] += AA
        topScores.append(temp)
    
    # For each candidate string, normalize the candidate by taking the average of each amino acid + the 
    # corresponding data point. 
    # Eliminate candidates that have 3 or fewer candidate amino acids. The assumption is that the data
    #  point will have at least three true amino acid matches, with a 4th one included under consideration
    #  to further eliminate candidates that arise from random chance.
    #  Amino acid weights will be used in later analysis. If fewer than 3 true amino acids appear in the data,
    #  it is likely the program will fail to predict the cyclopeptide.
    finalGroups = []
    for seq in topScores:
        matchingDataPoints = []
        for AA in seq[1]:
            matchingDataPoints.append( [ ApproxSet( seq[0]-AminoAcids_String_Weight[AA], set(esData), tolerance = 1.000000001 ), AA ] )
        finalScore = dict()
        for dp in matchingDataPoints:
            m = ApproxSet( dp[0]+AminoAcids_String_Weight[ dp[1] ], set(finalScore.keys()) )
            if m == -1:
                finalScore[dp[0]+AminoAcids_String_Weight[ dp[1] ]] =  dp[1] 
            else:
                tempM = m
                tempM *= len(finalScore[m])
                tempM += dp[0]+AminoAcids_String_Weight[ dp[1] ]
                tempO = finalScore[m] + dp[1] 
                tempM /= float( len( tempO ) )
                finalScore.pop(m)
                finalScore[tempM] = tempO
        for i in finalScore:
            if len(finalScore[i]) > 3:
                finalGroups.append( [ i, finalScore[i] ] )
    return finalGroups 



##################################
#
#
#   Section 4:
#   Sequence Cyclopeptide Part 1
#   The following methodology
#   is effective at sequencing
#   MOST of the cyclopeptide
#
#
##################################

##################
# Some background for this calculation. As covered in the CreateConnectivityGraph() comments, we will have a series of 
#  couplet values that could represent another junction in the cyclopeptide. Each true coupling
#  represents a unique junction compared to the other couplets. When taken together, one could
#  theoretically reconstruct the linear peptide that represents the segemnts joining each of these
#  junctions.
# NOTE: There are likely several ways to do this. The route I chose was the one that seemed to work
#  best out of the ones I tried, but we have officially entered speculative, adaptable territory. Each choice
#  leads to another choice which leads to another choice, at any of which multiple options are technically "right,"
#  though some may work better than others. In essence, here are the puzzle pieces, how they come together is up
#  to you.
#################

# This would be easier to explain visually, but I'll do my best here.
#  Each true data point represents 2 junctions of the original cyclopeptide.
#  You can use this to calculate the liklihood that any other two data points
#  conjoin this data point at the edges.
#
# dp is the data point value being tested (dp1 in the following explanation). This should be
#  a smaller data point to increase the liklihood over it being part of an overlap.
# esData is the experimental spectrum data.
# full is the predicted weight of the full cyclopeptide.
#
# Returns cGraph, a connectivity graph (see explanation)
# Returns flipAndAdd, a collection of couplets corresponding to single junctions (see explanation)
def CreateConnectivityGraph( dp, esData, full ): 
    # The "flip and add" technique is fairly straightforward. Take a data point (dp1).
    #  Flip any other data point (dp2) to it's opposite strand. Add the first data point to
    #  the opposite strand. If said strand also happens to match a third data point (dp3), you
    #  just found two data points (dp2 and dp3) that overlap on one side and meet somewhere else
    #  in the cyclopeptide (at a third junction). The overlap is the original data
    #  point (dp1) being tested. Thus, each "flip and add" couplet that corresponds to true data points
    #  represents two data points that meet at a third data point.
    sameSideSmall = CreateSmallSegments( esData, full )
    doubleData = DoubleTheData( sorted(sameSideSmall), full )
    sameSideLarge = []
    for i in sameSideSmall:
        sameSideLarge.append(full-i)
    flipAndAdd = []
    for i in sorted(sameSideSmall):
        flip = full - i
        add = flip + dp
        for j in sameSideLarge:
            if Approx( add, j ):
                flipAndAdd.append( [round(i,2), round(j,2)] )

    # Each true "flip and add" couplet represents a third junction. Furthermore, both "flip and add"
    #  couplet values are derived from the same two junctions (those represented by dp1). So, between
    #  each true "flip and add" couplet, every value pairs with one of the two from another couplet in sharing
    #  a junction from dp1. Thus, subtracting one value of one couplet from one value of another couplet could represent
    #  a third datapoint. Furthermore, each couplet's opposite would achieve the same effect.
    # The connectivity graph represents both values achieved from subtracting a value from one couplet from
    #  both values in another couplet. Both x and y axis are the length of the flip and add values.
    cGraph = [ [[0] for i in xrange(len( flipAndAdd ))  ] for i in xrange( len( flipAndAdd ) ) ]
    for i in xrange(len(flipAndAdd)):
        for j in xrange(len(flipAndAdd)):
            segments = []
            seg1 = abs(flipAndAdd[i][0]-flipAndAdd[j][0])
            if seg1 > 50:
                segments.append( round(seg1,2) )
            seg2 = abs(flipAndAdd[i][0]-flipAndAdd[j][1])
            if seg2 > 50:
                segments.append( round(seg2,2) )
            if len(segments) != 0:            
                cGraph[i][j] = segments

    return cGraph, flipAndAdd    

# Generates k4mers from flipAndAdd couplet values. 
#  The k4mers are calculated in a particular format. As indicated by the name, there are 4 segments.
#  The first segment is the original data point(dp1). The second segment is the smaller section minus the
#  original data point. The third segment represents the outer edges of the linear peptide we are
#  trying to reconstruct. The fourth represents what's left to form the cyclopeptide. In that sense, there are only
#  two of these k4mers that accurately represent the complete cyclopeptide.
#   NOTE: You may want to calculate the new segments by taking the reverse
#   of a segment as opposed to subtracting the dp. I suspect it will be more
#   accurate that way.
#
# Takes the flipAndAdd couplet values
# Takes the original data point (dp1)
# Takes the predicted full weight of the cyclopeptide.
#
# Returns a list of tuples. The first element of each tuple is the k4mer. The second is the 
#  flipAndAdd indicies that were used to create that k4mer, indicies that match the connectivity Graph
#  for later analysis.
def GenerateLength4Kmers( FA, dp, full ):
    k4mers = []
    for i in xrange( len(FA) ):
        for j in xrange( i+1, len(FA) ):
            seg3 = FA[j][0]-FA[i][0]
            seg2 = FA[i][0] - dp
            seg4 = FA[j][1] - dp
            seg3 = round(seg3,2)
            seg2 = round(seg2,2)
            seg4 = round(seg4,2)
            if seg2 > 50 and seg3 > 50 and seg4 > 50:
                k4mers.append( [ [dp, seg2, seg3, seg4], [i,j] ] )
    return k4mers

# Serves two functions. It eliminates chains who have non-amino acid
#  segments who's weight is under 200. Furthermore, for segments that
#  Are approximately the value of an amino acid, it changes the segment
#  to match the amino acid weight exactly, thereby normalizing the 
#  chain segment.
#
# Takes list of cyclopeptide chains.
#
# Returns Trimeed and normalized list of cyclopeptide chains.
def TrimByAAs( kmers ):
    for i in xrange( len(kmers)-1,-1,-1 ):
        deleteThem = False
        for j in xrange(len(kmers[i][0])):
            match = ApproxSet( kmers[i][0][j], set(AminoAcids_String_Weight.values()) )
            if kmers[i][0][j] < 200 and match == -1:
                deleteThem = True
            elif kmers[i][0][j] < 200 and match != -1:
                kmers[i][0][j] = match
        if deleteThem:
            kmers.pop(i)
    return kmers

# Runs the Trim algorithm on a different leaderboard format (keeps flipAndAdd indicies paired
#  with leaderboard chain)
def TrimWithPoints( leaderboard, esData, N ):
    if N > len(leaderboard):
        return leaderboard
    scores = [ 0 for i in leaderboard ]
    for i in xrange( len(leaderboard) ):
        scores[i] = Score( leaderboard[i][0], set(esData) )
    scores, leaderboard = (list(x) for x in zip(*sorted(zip(scores, leaderboard), key=lambda pair: pair[0])))
    for i in xrange( len(scores)-1,-1,-1 ):
        if scores[i] < scores[len(scores)-N]:
            leaderboard = leaderboard[i+1:]
            break
    return leaderboard

# If any of the new kmers have a higher score than the previous leaders,
#  the new highest-scoring kmers replace the previous leaders.
#
# takes kmers, a list of kmers 
# takes esData, the experimental spectrum
# takes previousLeaders, the leaders prior to this function
#
# Returns the list of new leaders
def FindLeaderPeptide( kmers, esData, previousLeaders ):
    newLeaders = previousLeaders
    highScore = Score( newLeaders[0], esData )
    for i in kmers:
        if Score( i, esData ) > highScore:
            newLeaders = [i]
            highScore = Score( i, esData )
        elif Score( i, esData ) == highScore:
            newLeaders.append(i)
    return newLeaders

# Checks the linear chain against values in graph.
#  Every segment of the linear peptide we are trying to find is located in the connectivity graph.
#  Thus, if any segment does not match a value in Graph (as corresponded by the order of indicies given)
#  it is considered invalid.
# 
# Takes chain (a linear chain)
# Takes points, the list of indices represented by chain
# Takes Graph, the connectivity graph calculated earlier
#
# Returns True if all elements of the linear theoretical spectrum are accounted for, False if otherwise.
def CheckIfChainFitsAllData( chain, points, Graph ):
    tsc = LinearTheoreticalSpectrum( chain )
    for i in xrange( len(points)-1 ):
        for j in xrange(i+1, len(points) ):
            values = Graph[points[i]][points[j]]
            for k in values:
                if Approx( tsc[0],k, tolerance = 0.600000001 ):
                    tsc.pop(0)
                    break
    if len(tsc) == 0:
        return True
    else:
        return False

# Extends a list of kmers by inserting a new possible fragment.
#  To put in broader terms, I'm trying to add the "next" segment of the linear
#  peptide represented by the values in flipAndAdd. Thus, the minimum value
#  of each graph cell is added, excluding indices that correspond to flipAndAdd
#  values that have already been included. 
#
# takes kmers, a list of kmers (in the [kmer, FA indices] format)
# Takes Graph, the connectivity graph
#
# returns a list of extended kmers (including updated FA indices)
def ExtendKmers( kmers, Graph ):
    newKmers = []
    for kmer in kmers:
        for nextPoint in xrange( len(Graph) ):
            if nextPoint not in kmer[1]:
                #necessary????
                numAddedSegments = len(kmer[1])-1
                tempLinearChain = kmer[0][2:]
                newSegment = min( Graph[ kmer[1][-1] ][nextPoint])
                endSegment = tempLinearChain[-1]-newSegment
                endSegment = [round(endSegment, 2)]
                tempLinearChain = tempLinearChain[:-1] + [newSegment]  
                newPoints = [i for i in kmer[1]]
                newPoints.append(nextPoint)
                if CheckIfChainFitsAllData( tempLinearChain, newPoints, Graph ):
                    newChain = kmer[0][0:2] + tempLinearChain + endSegment
                    newKmers.append( [ newChain, newPoints ] )
    return newKmers

# Tries to recreate the cyclopeptide represented by the original data point (dp1).
#  Segments are added until all other options either result in a non-amino acid segment
#  that is less than 200 appears or all FA values have been exhausted (the former is much
#  more likely than the latter).
# The resulting may not be 100% accurate representations of the original cyclopeptide, but 
#  they will be close. If you generate several cyclopeptides that are "almost" accurate,
#  and that differ in different ways, a consensus chain could be created where overlapping
#  values are discovered between candidate cyclopeptides.That will come later, though.
#
# Takes kmers, list of kmers (specifically k4mers)
# Takes esData, the experimental spectrum
# Takes Graph, the connectivity graph
#
# Returns cyclopeptide candidate(s) that score the highest.
def GenerateLeaderPeptides( kmers, esData, Graph ):
    leaders = [[0,0]]
    kmers = TrimByAAs( kmers )
    while len(kmers)!= 0:
        leaders = FindLeaderPeptide( [i[0] for i in kmers], esData, leaders )
        kmers = ExtendKmers( kmers, Graph )
        kmers = TrimByAAs( kmers )
    leaders = EliminateDuplicateChains( leaders )
    return leaders

# Scores the alignment of two chains. Each shared junction between the two chains (including the first one)
#  results in a higher score.
#
# Takes chain1, a cyclopeptide chain.
# Takes chain2, a cyclopeptide chain.
#
# Returns the alignment score of the two cyclopeptides.
def ChainsAlignmentScore( Chain1, Chain2 ):
    score = 0
    chain1 = [i for i in Chain1]
    chain2 = [i for i in Chain2]
    while len(chain1) != 0 and len(chain2) != 0:
        if Approx(sum(chain1),sum(chain2)):
            score += 1
            chain1.pop(0)
            chain2.pop(0)
        elif sum(chain1) > sum(chain2):
            chain1.pop(0)
        elif sum(chain2) > sum(chain1):
            chain2.pop(0)
    return score 

# Takes two links and finds the best possible alignments between them.
#
# Takes chain1, a cyclopeptide chain
# Takes chain2, a cyclopeptide chain
#
# Returns two objects. A list of best alignments to chain1 (alternates of chain2)
#  and the score of those alignments. THIS MEANS that the best alignments
#  DO NOT contain any variations of chain1.
def AlignChains( chain1, chain2):
    bestAlignments = []
    bestAlignmentScore = 0
    rChain2 = chain2[::-1]
    for i in xrange( len(chain2) ):
        x = chain2[i:] + chain2[0:i]
        y = rChain2[i:] + rChain2[0:i]
        score1 = ChainsAlignmentScore( chain1, x )
        if score1 > bestAlignmentScore:
            bestAlignmentScore = score1
            bestAlignments = [x]
        elif score1 == bestAlignmentScore:
            bestAlignments.append(x)
        score2 = ChainsAlignmentScore( chain1, y )
        if score2 > bestAlignmentScore:
            
            bestAlignmentScore = score2
            bestAlignments = [y]
        elif score2 == bestAlignmentScore:
            bestAlignments.append(y)
    bestAlignments = EliminateDuplicateChains( bestAlignments )
    return bestAlignments,bestAlignmentScore 

# Takes two chains and merges them conservatively. By conservatively,
#  I mean that the final product only includes fragments that match
#  both chains. This, chain1 could have 6 segments, and chain2 could have
#  4, but if only three junctions are shared then the finished product
#  will have 3 segments. The chains MUST BE ALIGNED IN ADVANCE.
#
# Takes Chain1, a cyclopeptide chain
# Takes Chain2, another cyclopeptide chain
#
# Returns the merged final chain.
def MergeMatchedChainsConservative( Chain1, Chain2 ):
    finalChain = []
    chain1 = [i for i in Chain1]
    chain2 = [i for i in Chain2]
    while len(chain1)!= 0 and len(chain2)!= 0:
        if Approx( chain1[0],chain2[0] ):
            finalChain.append(round( (chain1[0]+chain2[0])/2,2 ))
            chain1.pop(0)
            chain2.pop(0)
        elif chain1[0] > chain2[0]:
            if len(chain2) > 1:
                chain2[1] += chain2[0]
            chain2.pop(0)
        else:
            if len(chain1) > 1:
                chain1[1] += chain1[0]
            chain1.pop(0)
    return finalChain

# Applies MergeMatchedChainsConservative() to the elements in a list.
#  chain1 is used as the standard alignment.
# 
# Takes chain1, a cyclopeptide sequence
# Takes otherChains a list of cyclopeptide sequences that have been aligned to chain1
#
# Returns a list of merged chains.
def CreateConsensusSequenceSmallScale( chain1, otherChains ):
    finalChains = []
    for i in otherChains:
        chain = MergeMatchedChainsConservative( chain1, i )
        finalChains.append(chain)
    return finalChains

# This function comprehensively uses all of the preceding functions in the section, and
#  after CreateConnectivityGraph() consider it to be the second main step in predicting
#  the cyclopeptide sequence.
# In summary, one can generate at least one cyclopeptide that mostly lines up with the original cyclopeptide
#  from a single amino acid. This function calculates said cyclopeptide given a particular amino acid.
#
# Takes AminoAcid, a single-character string that has a matching value in the AminoAcid_String_Weight dictionary
# Takes full, the predicted full weight of the cyclopeptide
# Takes esData, the experimental spectrum data set.
#
# Returns any consensus chain(s). These chains each could represent some portion of the cyclopeptide.
#  It's possible that multiple consensus chains will be generated, and that the different chains hold different
#  aspects of the original cyclopeptide. Note: the summary weight of each chain should be approximately the 
#  weight of the full cyclopeptide. However, the chain may have 6 segments, 5 of which match up with the
#  full cyclopeptide. That kind of thing.
def GenerateConsensusSequencesFromAminoAcid( AminoAcid, esData, full ):
    g,FA = CreateConnectivityGraph( AminoAcids_String_Weight[AminoAcid], esData, full )
    
    k4mers = GenerateLength4Kmers( FA, AminoAcids_String_Weight[AminoAcid], full )
    leaderPeptides = GenerateLeaderPeptides( k4mers, esData, g )
    alignedChains = []
    consensus = []
#NOTE: this may be necessary in the future
#    if len(leaderPeptides) == 1:
#        consensus = leaderPeptides
    for i in xrange( len(leaderPeptides) ):
        for j in xrange( i+1, len(leaderPeptides) ):
            alignedChains,s = AlignChains( leaderPeptides[i],leaderPeptides[j] )
            consensus += CreateConsensusSequenceSmallScale( leaderPeptides[i], alignedChains )
    consensus = EliminateDuplicateChains( consensus )
    consensus = [i for i in consensus if len(i) != 0]
    return consensus

# Similar to MergeMatchedChainsConservatively(), however, the end result is quite different. Presuming
#  the two chains are matched as well as can be, the resulting chain will have more segments than either
#  chain (presuming that one is not totally included in the other). For example, if 3 junctions are shared
#  between chain1 and chain2, but chain1 has 5 junctions and chain2 has 6, the resulting chain will have
#  8 (shared + unique1 + unique2).
#
# Takes chain1, a cyclopeptide chain.
# Takes chain2, a cyclopeptide chain that has been matched to chain1.
#
# Returns a fragmented cyclopeptide chain.
def MergeMatchedChainsDestructive( Chain1, Chain2 ):
    finalChain = []
    chain1 = [ i for i in Chain1 ]
    chain2 = [ i for i in Chain2 ]
    while len(chain1)!= 0 and len(chain2)!= 0:
        if len(chain1) == 0:
            finalChain += chain2
        elif len( chain2 ) == 0:
            finalChain += chain1
        else:
            if Approx( chain1[0],chain2[0] ):
                finalChain.append(round( (chain1[0]+chain2[0])/2,2 ))
                chain1.pop(0)
                chain2.pop(0)
            elif chain1[0] > chain2[0]:
                finalChain.append(chain2[0])
                chain1[0] -= chain2[0]
                chain2.pop(0)
            else:
                finalChain.append( chain1[0] )
                chain2[0] -= chain1[0]
                chain1.pop(0)
    return finalChain

# In a later function, the distances between possible amino acids will be calculated.
#  Each distance will be stored in a graph (aaGraph). As can be expected, not all distances
#  between possible amino acids are accurately portrayed in the incomplete cyclopeptides.
#  As such, each possible distance between possible amino acids is scored. The score is
#  based on whether or not the distance is verified by another amino acid distance.
#  The distances with the highest score are kept, all others are removed.
#
# Takes dp1, an index that represents an amino acid in aaGraph and possAminoAcids
# Takes dp2, a similar index but for a different amino acid
# aaGraph, a graph containing distances between amino acids
# possAminoAcids, a list containing amino acid weights (corresponds to indices)
#
# Returns the highest-scoring values in the graph cell represented by [dp1, dp2]
def TrimCellContents( dp1, dp2, full, aaGraph, possAminoAcids ):
    cell = aaGraph[dp1][dp2]
    score = [0 for i in cell]
    for i in xrange( len( aaGraph ) ):
        if i != dp1 and i!= dp2:
            dp2ValueSet = set()
            for x in aaGraph[dp2][i]:
                dp2ValueSet.add( x )
                xR = full-x-possAminoAcids[i]-possAminoAcids[dp2]
                dp2ValueSet.add(xR)
            for x in xrange( len( cell ) ):
                for y in aaGraph[dp1][i]:
                    temp = cell[x] + possAminoAcids[dp1] + y
                    if ApproxSet( temp, dp2ValueSet )!= -1:
                        score[x] += 1
    newCell = [i for i in cell]
    if max(score) == 0:
        return []
    for x in xrange( len(score)-1,-1,-1 ):
        if score[x] != max(score):
            newCell.pop(x)
    return newCell

# This is a function that generates all possible permutations of a list of 
#  indices with no repeats and omissions allowed. It'll be used later.
#
# Takes indices, a list of integers corresponding to indices
#
# Returns a list of lists containing all permutations of indices.
def GeneratePermutations( indices ):
    combos = [ [i] for i in indices ]
    overall = [ i for i in combos ] 
    final = 1
    while final < len(indices):
        temp = []
        final += 1
        for c in combos:
            for i in indices:
                if i not in c:
                    temp.append( c+[i] )
                    overall.append( c+[i] )
        combos = temp
    for i in xrange( len(overall) ):
        temp = overall[i]
        temp = sorted(temp)
        overall[i] = temp
    overall = sorted(overall)
    for i in xrange( len(overall)-1,-1,-1 ):
        if i != 0 and overall[i] == overall[i-1]:
            overall.pop(i)
    return overall

# Generates cyclopeptides based on a series of connections between amino acids.
#  Specifically, this algorithm will be run on each "candidate" amino acid. All other amino acids
#  That have a possible relation to this amino acid are included in the analysis. Because there are
#  Two possible orientations for each amino acid in relation to the candidate amino acid 
#  (clockwise or counterclockwise, neither of which we know), the algorithm builds the cyclopeptide
#  one amino acid at a time.
# It is distinctally possible that some of the possibleRelations do not actually relate to the candidate
#  amino acid. As such, this algorithm will be run on all possible combinations of possible relations.
#  That is the purpose of the GeneratePermutations() algorithm - all possible permutations of indices
#  that relate to possibleRelations are calculated, and this function is applied to each one.
#  Sequences are then scored and either discarded or kept.
#
# Takes indices, a result from the GenereatePermutations() function, which was run on possibleRelations indicies
# Takes possibleRelations, a list of triplet lists. [0] is the amino acid, [1] is the smaller distance between
#  the candidate amino acid and the [0] amino acid, and [2] is the opposite distance of [1]. if [1] is 0, that
#  indicates that the amino acids could be neighboring one another.
# Takes originalCyclopeptide, a two-part list. [0] is the candidate amino acid, [1] is the rest of the cyclopeptide.
# Takes aaGraph, a graph of amino acid distances.
# Takes aaIndex, a list of amino acid values that corresponds to aaGraph.
# Takes full, the full sequence.
# 
# Returns the list of all possible cyclopeptides of this particular permutation that were not cut. 
def GenerateNewCyclopeptides( indices, possibleRelations, originalCyclopeptide, aaGraph, aaIndex, full ):
    cyclos = [originalCyclopeptide]
    for i in indices:
        temp = []
        for c in cyclos:
            if possibleRelations[i][1] != 0:
                newCyclo1 = [ originalCyclopeptide[0], possibleRelations[i][1], possibleRelations[i][0], possibleRelations[i][2] ]
            else:
                newCyclo1 = [ originalCyclopeptide[0], possibleRelations[i][0], possibleRelations[i][2] ]
                
            newCyclo2 = [ originalCyclopeptide[0], possibleRelations[i][2], possibleRelations[i][0], possibleRelations[i][1] ]
            newCyclo1 = MergeMatchedChainsDestructive( c, newCyclo1 )
            newCyclo2 = MergeMatchedChainsDestructive( c, newCyclo2 )
            if ScoreSequence( newCyclo1, aaGraph, aaIndex, full )!= -1:
                temp.append( newCyclo1 )
            if ScoreSequence( newCyclo2, aaGraph, aaIndex, full )!= -1:
                temp.append( newCyclo2 )
        cyclos = temp
    
    return cyclos

# The Scoring algorithm for each sequence. In looking over how it is used, it seems that only the top half
#  is utilized in this Algorithm. The top half essentially just determines that each segment that is below
#  200 daltons corresponds to an amino acid, and if it doesn't the cyclopeptide is given a -1 score, which is a
#  signal that it should be cut. NOTE I'm not sure if the bottom half is utilized, if it's not I should cut it out
#  or find a use for it.
#
# Takes sequence, a cyclopeptide chain candidate
# Takes aaGraph, a graph of amino acid distances
# Takes aaIndex, a list of amino acid values that corresponds to aaGraph
# Takes full, the full sequence.
#
# Returns a score. If score is -1, the sequence is considered invalid.
def ScoreSequence( sequence, aaGraph, aaIndex, full ):
    indices = [ -1 for i in sequence ]
    score = 0
    for i in xrange( len(sequence) ):
        if sequence[i] < 200 and ApproxSet( sequence[i], set(AminoAcids_String_Weight.values() ) )==-1:
            return -1
        elif sequence[i] in aaIndex:
            indices[i] = aaIndex.index( sequence[i] )
    for i in xrange( len(indices) ):
        for j in xrange( i+1, len(indices) ):
            if indices[i] != -1 and indices[j] != -1:
                cellSet = set()
                for x in aaGraph[indices[i]][indices[j]]:
                    cellSet.add( x )
                    cellSet.add( full-x-aaIndex[indices[i]]-aaIndex[indices[j]] )
                inBetween = sum( sequence[i+1:j] )
                if ApproxSet( inBetween, cellSet )!= -1 or ( ApproxSet( inBetween, cellSet )== 0 and inBetween == 0 ):
                    score += 1
    
    return score

# Generates a cyclopeptide sequence candidate for a specific amino acid. The idea is to calculate
#  a cyclopeptide through the relation of each amino acid to the original "candidate" amino acid.
#  Many possible cyclopeptides are generated, as each amino acid relation could be in either of
#  two directions when including multiple amino acids. The ones with the highest score in relation
#  to the original data set are returned.
#
# Takes index, the index of the amino acid being analyzed.
# aaGraph, a graph of amino acid distances.
# aaIndex, a list of amino acid values that corresponds to aaGraph (and index)
# Takes full, the full sequence.
# Takes esData, the experimental spectrum data set.
# 
# Returns topCyclos, a list of the top cyclopeptides produced.
def SequenceCyclopeptide( index, aaGraph, aaIndex, full, esData ):
    possibleRelations = []
    for i in xrange( len( aaGraph[index] ) ):
        if len(aaGraph[index][i]) != 0:
            for j in xrange( len(aaGraph[index][i]) ):
                possibleRelations.append( [ aaIndex[i], aaGraph[index][i][j], round(full-aaGraph[index][i][j]-aaIndex[i]-aaIndex[index],2) ] )
    indices = [i for i in xrange( len( possibleRelations ) )]
    allPermutations = GeneratePermutations( indices )
    overall = []
    for i in allPermutations:
        overall += GenerateNewCyclopeptides( i, possibleRelations, [ aaIndex[index], full-aaIndex[index] ], aaGraph, aaIndex, full )
    overall = EliminateDuplicateChains(overall)
    
    for i in xrange(len(overall)-1,-1,-1):
        if ScoreSequence( overall[i], aaGraph, aaIndex, full ) == -1:
            overall.pop(i)
    topCylcos = []
    highScore = -1
            
    for i in overall:
        if Score( i, set(esData) ) > highScore:
            highScore = Score(i, set(esData))
            topCyclos = [i]
        elif Score( i, set(esData )) == highScore:
            topCyclos.append(i)
    return topCyclos

# In evaluating the various cyclopeptides that resulted from SequenceCyclopeptide() as used on each
#  possible amino acid, the real sequence appeared multiple times while all others did not. This,
#  this function scores each cyclopeptide by the number of times it gets repeated in the list, starting
#  from the back. Later, the most frequently repeated sequence will be chosen as the ideal cyclopeptide.
#
# Takes chains, a list of cyclopeptide candidates.
#
# Returns scores, a list of scoring integers. The index of a score relates to the index of the cyclopeptide.
def MostFrequentlyRepeatingChains( chains ):
    scores = [0 for i in chains]
    for i in xrange( len(chains)-1,-1,-1 ):
        for j in xrange( i-1,-1,-1 ):
            if len(chains[i])==len( chains[j] ):
                for k in xrange(len(chains[j])):
                    newChain = list(chains[j][k:]) + list(chains[j][0:k])
                    reverseChain = newChain[::-1]
                    if AreChainsSame( chains[i], newChain ) or AreChainsSame( chains[i], reverseChain ):
                        scores[i] += 1
    print scores
    return scores

# Runs everything described above, pops out the most frequently repeated cyclopeptide chain as the "ideal"
#  cyclopeptide.
#
# Takes esData, the experimental spectrum data set
# Takes aminoAcidString, a string of amino acids likely to be a part of the full cyclopeptide (see section 3)
# Takes full, the full sequence.
#
# Returns a candidate cyclopeptide.
def CreateConsensusSequenceLargeScale( esData, aminoAcidString, full ):
    consensus = []
    for i in xrange( len( aminoAcidString ) ):
        c = GenerateConsensusSequencesFromAminoAcid( aminoAcidString[i], esData, full )
##########
        print aminoAcidString[i]
        for i in c:
            print i
        consensus += c
    aaLinks = []
    for chain in consensus:
        for i in xrange( len(chain) ):
            if ApproxSet( chain[i], set(AminoAcids_String_Weight.values()) )!= -1:
                for j in xrange( i+1, len(chain) ):
                    if ApproxSet( chain[j], set( AminoAcids_String_Weight.values() )) != -1 :
                        diff = sum( chain[i+1:j] )
                        diff2 = full - (diff + chain[i] + chain[j])
                        if diff2 < diff:
                            diff = diff2
                        if Approx( diff, 0 ):
                            diff = 0
                        aaLinks.append( [ min( chain[i],chain[j] ), max( chain[i],chain[j] ), round(diff,2) ] )
#                        aaLinks.append( [ max( chain[i],chain[j] ), min( chain[i],chain[j] ), diff ] )
    possibleAAs = set()
    for i in aaLinks:
        for j in xrange(2):
            possibleAAs.add( i[j] )
    print "possibleAAs: " + str(sorted(possibleAAs))
    aaGraph = [ [ [] for i in possibleAAs ] for i in possibleAAs ]
    aaIndex = sorted(possibleAAs)
    for link in aaLinks:
        aaGraph[ aaIndex.index( link[0] ) ][ aaIndex.index( link[1] ) ] += [ link[2] ]
        aaGraph[ aaIndex.index( link[1] ) ][ aaIndex.index( link[0] ) ] += [ link[2] ]
    for i in xrange( len( aaGraph ) ):
        for j in xrange( len(aaGraph) ):
            newCell = [k[0] for k in EliminateApproximateDuplicateValues( aaGraph[i][j] )]
            aaGraph[i][j] = sorted(newCell)#EliminateApproximateDuplicateValues( aaGraph[i][j] ))
#########
#    print "\n"
#    for i in xrange( len(aaIndex)):
#        print str(i) + ": " + str(aaIndex[i])
#    print "\n"
    for i in xrange( len(aaGraph) ):
        for j in xrange( i+1, len( aaGraph ) ):
            if len( aaGraph[i][j] )>1:
                aaGraph[i][j] = TrimCellContents( i, j, full, aaGraph, aaIndex )
                aaGraph[j][i] = aaGraph[i][j]
    for i in aaGraph:
        print i
##########
#    for i in xrange( len(aaGraph)):
#        for j in xrange( i+1, len( aaGraph ) ):
#            if len( aaGraph[i][j] )!= 0:
#                print str(i) + " " + str(j)+ " " + str(aaGraph[i][j])
   
    topCyclos = []
    print "length of aaGraph: " + str(len(aaGraph))
    for i in xrange( len(aaGraph) ):
        print "Generating topCyclos: " + str(i)
        topCyclos += SequenceCyclopeptide( i, aaGraph, aaIndex, full, esData )
    
    print "top:"    
    mostF = MostFrequentlyRepeatingChains( topCyclos ) 
    print "TopCyclos for: " + str(full)
    for i in topCyclos:
        print i
    print "\n"
    if len(mostF) == 0:
        return []
    topCyclo = topCyclos[mostF.index(max(mostF))]
    if max(mostF) == 0:
        return []
    return topCyclo






##################################
#
#
#   Running The Algorithm
#
#
##################################



# For Section 3, predicting the full sequence of the cyclopeptide
#  In testing mode, will be more automated at the end
fullSizes= PredictFullSizes( data, 5, 100 )
for i in xrange( len(fullSizes) ):
    print str(i) + ": " + str(fullSizes[i]) 

seqIndex = raw_input("Which sequence would you like to run the cyclopeptide analysis on? ")
seqIndex = int(seqIndex)
#seqIndex = 0

originalFull = fullSizes[seqIndex][0]
originalAAs = fullSizes[ seqIndex ][1]

#originalFull = 1322.7
#originalAAs = "AFMVXZ"
#originalAAs = "AVXZMF"

print originalFull
print originalAAs

# For Section 4, partially sequencing the cylcopeptide 
rData = [ i for i in data if i < originalFull-50  ]
#for i in originalAAs:
#    rData.append( originalFull - AminoAcids_String_Weight[i] )

#c1 = GenerateConsensusSequencesFromAminoAcid( "X", rData, originalFull )
dd = DoubleTheData( rData, originalFull )

diff = -0.2
total = []
for i in xrange( 3 ):
    diff += 0.1
    c = CreateConsensusSequenceLargeScale( data, originalAAs, originalFull+diff )
    total.append(c)

for i in total:
    print i
#print "C: " + str(c)
