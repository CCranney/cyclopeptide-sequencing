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

tempCyclos = []
fullCyclos = []

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
#  If the data point is a true represenBtative of the cyclopeptide, then full-dp
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
#  to represeEnting the original cyclopeptide, the more it's own theoretical spectrum
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
#   5 parts of the cyclopeptide
#
#
##################################

# This algorithm subtracts every data point from every smaller
#  data point to see of the result matches a third (smallest)
#  data point. If there is a match, the data points are scored.
#  After ordering the data points by score, the top N data points
#  are selected. I've generally picked 10.
#
# Ultimately, the 10+ data points that the function outputs
#  generally have a higher liklihood of being part of the theoretical
#  spectrum of the full cyclopeptide. 
#
# esData is the experimental data set or whatever other data set
#  one wishes to run the algorithm on. This is in sorted list 
#  format.
# 
# Returns a list of 10+ data points.

def MakeTopN( esData, N ):
    connections = []
    score = dict( zip( esData, [ 0 for i in xrange( len( esData ) ) ] ) )
    for i in xrange( len(esData)-1, -1, -1 ):
        for j in xrange( i-1, -1, -1 ):
            for k in xrange( j-1, -1, -1 ):
                if esData[i]-esData[j] > esData[k]-0.3 and esData[i]-esData[j] < esData[k]+0.3:
                    connections.append( [ esData[i] , esData[j] , esData[k] ] )
                    score[esData[i]] += 1
                    score[esData[j]] += 1
                    score[esData[k]] += 1
    sorted_score = sorted(score.items(), key=operator.itemgetter(1),reverse=True)
    cutoff = sorted_score[N-1][1]
    top10 = sorted([ i[0] for i in sorted_score if i[1] >= cutoff ])
    top10 = [i for i in top10]
    top10 = sorted(set( top10 ))
    top10 = [round(i,1) for i in top10]
    possAAs = [ ApproxSet(i, set(AminoAcids_String_Weight.values())) for i in top10 if ApproxSet( i, set(AminoAcids_String_Weight.values()) ) != -1 ]
    top10 = [ i for i in top10 if ApproxSet( i, set(AminoAcids_String_Weight.values()) ) == -1 ]
    top10 += [ AminoAcids_String_Weight[i] for i in originalAAs ] + possAAs
    top10 = sorted([ i[0] for i in EliminateApproximateDuplicateValues(top10)])

    return top10

# This function is at the heart of my algorithm.
#
# Takes a (small) datapoint. Then it runs through every
#  data point larger than it. When the large point is  
#  "flipped" to it's opposite strand, the smaLller
#  data point is "added" to the result, and the final
#  result is also in the dataset, you just found two
#  data points that may overlap over the smaller data
#  point and touch at a third junction that does not
#  border the smaller data point.
#
# Takes dp, the smaller data point.
# Takes esData, the data set to be used.
# Takes full, the full size of the cyclopeptide.
#  I will likely edit this out in favor of using
#  originalFull as a global variable.
#
# Returns a list of "FA" pairs corresponding to
#  the smaller data point.

def FlipAndAdd( dp, esData, full ):
    flipAndAdd = []
    doubleData = DoubleTheData( sorted(esData), full )
    doubleData = sorted(doubleData)
    for i in xrange( len(doubleData) ):
        flip = full - doubleData[i]
        add = flip + dp
        for j in xrange( i+1, len(doubleData) ):
            if Approx( add, doubleData[j]):
                flipAndAdd.append( [round(doubleData[j],2),round(doubleData[i],2) ] )
    return flipAndAdd

# This function merely edits a cyclopeptide to be
#  as accurate as possible. This is done by removing
#  0s and altering values to approximately similar
#  data points in the data set and/or amino acid weights.
# 
# At the present time it also edits out cyclopeptides
#  with non-amino acid values less than 200, which I will
#  need to edit out later when applying the algorithm to
#  cyclopeptides that have amino acids we don't know of.
#
# Takes cyclo, a list of floats corresponding to fragments
#  of a cyclopeptide.
#
# returns an edited cyclo.

def CorrectCyclo( cyclo ):
    for c in xrange(len(cyclo)-1,-1,-1):
        aa = ApproxSet( cyclo[c], set(AminoAcids_String_Weight.values()), tolerance = 0.600000000001 )
        match = ApproxSet( cyclo[c], set(dData), tolerance = 0.600000000001 )
        if Approx(cyclo[c],0):
            cyclo.pop(c)
        elif cyclo[c] < 200 and aa == -1:
            return []
        elif cyclo[c] < 200 and aa != -1:
            cyclo[c] = aa
        elif match != -1:
            cyclo[c] = match
    return cyclo

# This function first takes all the values generated
#  in the FA pairing for a particular small data point
#  (called "value"). These FA values may relate to the
#  smaller data point in that the smaller data point
#  is on one of the two ends of the FA values if they
#  represent the real cyclopeptide. So if you three
#  of the FA values involved overlap with one another
#  (As in, the differences between them are also in dData,
#  or and if you subtract the small data point it equals
#  another data point in the data set (using rData at the
#  present time)) then you can use those 3 values and the
#  smaller dataset to create a 5-part portion of the
#  cyclopeptide. 
#
# Takes value, a small and likely true data point
#
# Returns a list of 5mers

def Generate5mers( value ):
    k5mers = []
    fa = FlipAndAdd( value, dData, originalFull )
    fData = [i[0] for i in fa] + [ i[1] for i in fa ]
    fData = sorted(fData)

    for i in xrange( len(fData)-1,-1,-1 ):
        for j in xrange( i-1,-1,-1 ):
            for k in xrange( j-1,-1,-1 ):
                m = ApproxSet( fData[i]-value, set(rData) )
                n = ApproxSet( fData[j]-value, set(rData) )
                o = ApproxSet( fData[k]-value, set(rData) )
                p = ApproxSet( originalFull-value, set(rData) )
                check = True
                if ApproxSet(m-n, set(dData)) == -1 or ApproxSet( n-o, set(dData) )== -1:
                    check = False
                if m != -1 and n != -1 and o != -1 and p != -1 and check:
                    k5mers.append( [ round(value,1), round(originalFull-fData[i],1), round(fData[i]-fData[j],1), round(fData[j]-fData[k],1), round(fData[k]-value,1) ] )
    k5mers = [CorrectCyclo(i) for i in k5mers]
    k5mers = sorted( [ [Score(i,rData), i] for i in k5mers ], reverse = True )
    return k5mers

##################################
#
#
#   Section 5:
#   Sequencing Cyclopeptide Part 2 
#
#
##################################

# Takes a value and adds it to a cyclopeptide.
#  In this calculation, it is assumed that the
#  value (link) and the cyclopeptide are alreAady
#  matched. For example, if link = 5 and cyclo =
#  [ 1, 2, 7 ], the result would be [ 1, 2, 2, 5 ]
#
# Takes link, a float value
# Takes a cyclopeptide chain that the link
#  will be added to
#
# returns the new cyclopeptide if it does not
#  contain non-amino acids less than 200. If
#  it does, an empty list is returned.

def AddLinkToCyclo( link, cyclo ):
    for i in xrange( len(cyclo) ):
        temp = sum(cyclo[0:i+1])
        if Approx(link, temp):
            return cyclo
        if temp > link:
            part1 = link - (temp-cyclo[i])
            part2 = temp - link
            if part1 < 200:
                m = ApproxSet( part1, set( AminoAcids_String_Weight.values() ), tolerance = 0.6000001 )
                if m == -1:
                    return []
                part1 = m
            if part2 < 200:
                m = ApproxSet( part2, set( AminoAcids_String_Weight.values() ), tolerance = 0.6000001 )
                if m == -1:
                    return []
                part2 = m
            return CorrectCyclo(cyclo[0:i] + [part1, part2] + cyclo[i+1:])
    return [] 

# This algorithm relies on the assumption that, when
#  two small data points neighbor one another, any
#  FA pairings either data point creates that meet at
#  the same junction will have one value that is the
#  opposite fragment of a data point in the other
#  FA pairing. One can use that to create new links
#  in the cChain. This algorithm finds all possible
#  such junctions and returns the resulting kmers.
#
# Takes k5mer, an incomplete cyclopeptide representative.
#
# Returns cyclopeptide possibilities that are one length
#  longer than the input.

def AddLinkFromNeighbors( k5mer ):
    summary = set([ sum(k5mer[ 0:i ] ) for i in xrange( 1,len(k5mer) ) ] + [ originalFull ] + [0])
    FAs = [ FlipAndAdd( i, dData, originalFull ) for i in k5mer ]
    lengths = []
    for i in xrange(len(k5mer)):
        add = sum(k5mer[0:i])
        tempCyclo = k5mer[i:]+k5mer[0:i]
        neighbor1 = tempCyclo[0]
        neighbor2 = tempCyclo[-1]
        j = i-1
        if j < 0:
            j = len(k5mer)-1
        for k in xrange( len(FAs[j]) ):
            backwards = set( [ originalFull-FAs[j][k][0],originalFull-FAs[j][k][1] ] )
            for l in xrange( len(FAs[i]) ):
                m = ApproxSet( FAs[i][l][0], backwards, tolerance = 0.60000001 )
                n = ApproxSet( FAs[i][l][1], backwards, tolerance = 0.60000001 )
                if m != -1:
                    length = FAs[i][l][0] + add
                    if length > originalFull:
                        length -= originalFull
                    lengths.append(length) 
                if n != -1:
                    length = FAs[i][l][1] + add
                    if length > originalFull:
                        length -= originalFull
                    lengths.append(length) 

    lengths = EliminateApproximateDuplicateValues( lengths )
    for i in xrange( len(lengths)-1,-1,-1 ):
        if ApproxSet( lengths[i][0], summary, tolerance = 0.60000001 ) != -1:
            lengths.pop(i)
    lengths = [ AddLinkToCyclo( i[0], k5mer ) for i in lengths ]
    lengths = [ CorrectCyclo( i ) for i in lengths ]
    lengths = [ i for i in lengths if i != [] ]
    lengths = EliminateDuplicateChains( lengths )
    if len(lengths) == 1:
        return []
    return lengths 

# Determines if a cyclopeptide has all values that match an
#  amino acid, signifying that it is a complete cyclopeptide.
#
# takes cyclo, a length of floats
#
# returns True if all values correspond to amino acids,
#  False otherwise

def IsFullCyclopeptide( cyclo ):
    for peptide in cyclo:
        if ApproxSet( peptide, set( AminoAcids_String_Weight.values() ) )==-1:
            return False
    return True

# This algorithm needs work. Basically adds length after
#  length to the 5mers generated in the last part until
#  finally returning values that correspond to full
#  cyclopeptides.
#
# takes kmers, a list of partially complete cyclopeptides
#
# returns Nothing. Values are added to a global list called
#  tempCyclos.

def ConstructCyclopeptideRecursion( kmers ):
    if len(tempCyclos) >= 30:
        return None
    if len(kmers) == 0:
        return None 
    for kmer in kmers:
        if IsFullCyclopeptide(kmer):
            tempCyclos.append(kmer)
        else:
            newKmers = AddLinkFromNeighbors(kmer)
            ConstructCyclopeptideRecursion( newKmers )
    return None

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

originalFull = fullSizes[seqIndex][0] 
originalAAs = fullSizes[ seqIndex ][1]

print originalFull
print originalAAs

# For Section 4 and 5, sequencing the cyclopeptide
rData = [ i for i in data ]
for i in xrange(len(rData)-1,-1,-1):
    temp = originalFull-rData[i]
    if temp < 200 and ApproxSet( temp, set( AminoAcids_String_Weight.values() ), tolerance = 0.60000001 ) == -1:
        rData.pop(i)

dData = DoubleTheData( rData, originalFull )

top10 = MakeTopN(sorted(dData), 10)

# For testing purposes, since I keep overloading my recursion limit
top10 = [ AminoAcids_String_Weight[i] for i in originalAAs ]

for dp in sorted(top10):
    print dp
    tempCyclos = []
    t = Generate5mers( dp )
    ConstructCyclopeptideRecursion( [i[1] for i in t] )
    fullCyclos += tempCyclos
fullCyclos = EliminateDuplicateChains(fullCyclos)
fullCyclos = [ [ len(fullCyclos[i]), Score(fullCyclos[i],rData), fullCyclos[i] ] for i in xrange( len( fullCyclos ) ) ]
for i in sorted(fullCyclos):
    print i

