#!/usr/bin/python
# Victoria O'Neill and Matthew Farr
# CS 423, Lab 5
# fall 2015

######################################################################
# Determine the score of the optimal local alignment of two strings
# students: complete this function
######################################################################
def localAlignmentScore(s1, s2):

     # Scoring system
     MATCH = 5
     MISMATCH = -4
     GAP = -6

     # set table size
     NUM_ROWS = len(s2)+1
     NUM_COLS = len(s1)+1

     # Create table and fill it with zeros
     costs = createTable(NUM_ROWS, NUM_COLS, 0)

     # Create table for getting back the optimal alignment, fill table with "A"
     directions = createTable(NUM_ROWS, NUM_COLS, "A")
	
     # Fill the top row of the costs table with 0's, fill the top row of the
     # directions table with F's
     for i in range (0,NUM_COLS):
          costs[0][i] = 0
          directions[0][i] = "F"

     # Fill the left column of the costs table with 0's, fill the left column
     # of the directions table with F's
     for j in range (1, NUM_ROWS):
          costs[j][0] = 0
          directions[j][0] = "F"

     maxValue = 0
     maxRowPos = 0
     maxColPos = 0

     # Fill in each of the rows with alignment costs, starting at the top row
     # going left to right
     for y in range (1,NUM_ROWS):
          for x in range (1,NUM_COLS):
               # Calculate costs of a gap for the top and bottom sequences
               valTop = costs[y-1][x] + GAP
               valLeft = costs[y][x-1] + GAP

               # Increase score for a match, or decrease for a mismatch
               if (s1[x-1] == s2[y-1]):
                    valDag = costs[y-1][x-1] + MATCH
               else: valDag = costs[y-1][x-1] + MISMATCH

               # Calculate the maximum value between the three calculated
               # values and zero, and set that one as the cost for this
               # index of the table.
               val = max(valTop,valLeft,valDag,0)
               if (val > maxValue):
                    maxValue = val
                    maxRowPos = y
                    maxColPos = x
               costs[y][x] = val

               # Set the direction indicator table based on which value
               # was the maximum
               if val == 0:
                    # If no cost was greater than zero, set this as the
                    # beginning of a new substring
                    directions[y][x] = "F"
               # Otherwise, set the direction for where the maximum
               # value came from
               elif val == valLeft:
                    directions[y][x] = "L"
               elif val == valDag:
                    directions[y][x] = "D"
               else:
                    directions[y][x] = "T"
     
     # Print out table (only useful for small tables - used for debugging)
     # Comment out when you are satisfied that the algorithm is working
     printTable(costs, "costs.txt")
     printTable(directions, "directions.txt")

     # find optimal local alignment alignment
     align(directions, s1, s2, maxRowPos, maxColPos, "alignment.txt")

     # return optimal local alignment score
     return costs[maxRowPos][maxColPos]

####################################################################
# Create a 2D table with the given number of rows and columns
# and fills all entries with value given as a parameter
# (function completed for you)
####################################################################
def createTable(numRows, numCols, value):
     table = []
     row = 0
     # create 2D table initialized with value
     while (row < numRows):
          table.append([])    
          col = 0
          while (col < numCols):
               table[row].append(value)
               col = col + 1
          row = row + 1
     return table

##################################################################
# Print 2D table to file (only useful for small tables for short
# strings)
# Should have tabs between the values on each row
# Useful function for debugging purposes
# Complete this function
##################################################################
def printTable(table, filename):

     file = open(filename, 'w')

     # Get the size of the table
     row = len(table)
     col = len(table[0])

     # Print the table to a file, separating columns with tab characters
     # and separating rows with newline characters
     for y in range (0,row):
          for x in range (0,col):
               file.write(str(table[y][x]) + "\t")
          file.write("\n")
     file.close()

     return
	

################################################################
# Reconstruct the optimal alignment and print the alignment
# to a file. Because the sequences can be long, print the 
# alignment 50 characters on one line, the other string of 50 characters
# on the next line, and then skip one line, as follows:
# AATT--GGCTATGCT--C-G-TTACGCA-TTACT-AA-TCCGGTC-AGGC
# AAATATGG---TGCTGGCTGCTT---CAGTTA-TGAACTCC---CCAGGC
#
# TATGGGTGCTATGCTCG--T--TACG-CA
# TCAT--TGG---TGCTGGCTGCTT--ACA
#
# Complete this function
# direction is a 2D table, seq1 and seq2 are the original DNA
# sequences to align, and filename is the name of the output file
###############################################################
def align(direction, s1, s2, row, col, filename):

     file = open(filename, 'w')

     # Set x and y coordinates for traversing the direction table
     # at the the coordinates of the highest score in the costs
     # table. Set the alignment strings as empty strings
     x = col
     y = row
     currentDir = direction[y][x]
     topSeq = ""
     botSeq = ""

     # Traverse the table until the beginning of the substring
     # is encountered
     while currentDir != "F":

          # If the direction index points up, move up the table,
          # and align a nuceleotide from sequence 2 with a gap
          if direction[y][x] == "T":
               topSeq = "-" + topSeq
               botSeq =  s2[y-1] + botSeq
               y -= 1
          # If the direction index points left, move left in the table,
          # and align a nucleotide from sequence 1 with a gap
          elif direction[y][x] == "L":
               topSeq = s1[x-1] + topSeq
               botSeq = "-" + botSeq
               x -= 1
          # If the direction index points diagonally, move up and left
          # in the table, and align nucleotides from sequence 1 and 2
          # with each other
          elif direction[y][x] == "D":
               topSeq = s1[x-1] + topSeq
               botSeq = s2[y-1] + botSeq
               x-=1
               y-=1

          # Get next direction index from current position in the table
          currentDir = direction[y][x]

     # Print the alignment strings to the output file, in paired lines
     # of up to 50 characters each
     for i in range (0, len(topSeq), 50):
          file.write(topSeq[i:i+50] + "\n")
          file.write(botSeq[i:i+50] + "\n")
          file.write("\n")

     file.close()
	
     return


### End of functions ###################################


###################################################
### Testing #######################################
###################################################

# Calculate local alignment score of two sequences
#s = "AGCGTCTA"
#t = "TGCATCTCG"
#s = "GATTACA"
#t = "AATGACAC"
s = "AAGGTATGAATC"
t = "CAGTTGCAA"

optimalScore = localAlignmentScore(s, t)

print(s)
print(t)
print("Local alignment score: " + str(optimalScore))



