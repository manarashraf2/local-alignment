class ScoreMatrix:
    def __init__(self, match, mismatch,gap):
       self.gap = gap
       self.match = match
       self.mismatch = mismatch

def ZeroMatrix(row, column):
    return [[0 for row in range(len(row)+1)] for column in range(len(column)+1)]

def LocalAlignmentforDNA(sequence1, sequence2,score=ScoreMatrix(1,-2,-1)):

    matrix =ZeroMatrix(sequence1,sequence2)
    tmatrix = ZeroMatrix(sequence1,sequence2)
    maxScore = 0
    Index = (None, None)
    for i in range(1, len(sequence2) + 1):
        for j in range(1, len(sequence1) + 1):

            matchScore = score.match if sequence2[i - 1] == sequence1[j - 1] else score.mismatch
            diagonal = matrix[i - 1][j - 1] + matchScore
            up = matrix[i - 1][j] + score.gap
            left = matrix[i][j - 1] + score.gap
            matrix[i][j] = max(0, diagonal, up, left)
            if matrix[i][j] == 0:
                tmatrix[i][j] = 0
            elif matrix[i][j] == left:
                tmatrix[i][j] = 1
            elif matrix[i][j] == up:
                tmatrix[i][j] = 2
            elif matrix[i][j] == diagonal:
                tmatrix[i][j] = 3

            if matrix[i][j] >= maxScore:
                Index = (i, j)
                maxScore = matrix[i][j]
        for i in range(0, len(sequence2) + 1):
            print(matrix[i])
        print('\n')

    firstalignedseq = ""
    secondalignedseq = ""
    (I, J) = Index
    while tmatrix[I][J] !=0:
        if tmatrix[I][J] == 3:
            secondalignedseq += sequence2[I - 1]
            firstalignedseq += sequence1[J - 1]
            I -= 1
            J -= 1
        elif tmatrix[I][J] == 2:
            secondalignedseq += sequence2[I - 1]
            firstalignedseq += '-'
            I -= 1

        else:
            secondalignedseq += '-'
            firstalignedseq += sequence1[J - 1]
            J -= 1


    firstalignedseq = firstalignedseq[::-1]
    secondalignedseq = secondalignedseq[::-1]

    print(firstalignedseq)
    print(secondalignedseq)

def blosummatrix(matrix):
    file_lines = matrix.readlines()
    matrix.close()
    dictinary = {}
    protein_string = file_lines[0]
    protein_string = protein_string.split()

    i = 1
    while i <= (len(file_lines) - 1):
        rows = file_lines[i]
        rows = rows.split()

        j = 1
        for character in rows[1:25]:
            dictinary[protein_string[i - 1], protein_string[
                j - 1]] = character  # i,j changes (row, column) to (aa at row, aa at column) for keys
            j += 1
        i += 1

    return (dictinary)
def LocalAlignmentforProtein(x, y, Bdictinary):

    maxScore = 0
    Index = (None,None)
    gap = -1
    matrix = ZeroMatrix(x, y)
    tmatrix = ZeroMatrix(x, y)
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            if Bdictinary.get((x[j - 1], y[i - 1])) is None:
                score = int(Bdictinary.get((y[i - 1], x[j - 1])))
            else:
                score = int(Bdictinary.get((x[j - 1], y[i - 1])))
            #similarityScore = (matrix[i][0], matrix[0][j])
            up=int(matrix[i - 1][j]) + gap
            left=int(matrix[i][j - 1]) + gap
            diagonal = (matrix[i - 1][j - 1]) + score
            matrix[i][j]= max(0, diagonal, up, left)
            if matrix[i][j] == 0:
                tmatrix[i][j] = 0
            elif matrix[i][j] == left:
                tmatrix[i][j] = 1
            elif matrix[i][j] == up:
                tmatrix[i][j] = 2
            elif matrix[i][j] == diagonal:
                tmatrix[i][j] = 3
            if matrix[i][j] >= maxScore:
                Index = (i, j)
                maxScore = matrix[i][j]

        for i in range(0, len(y) + 1):
            print(matrix[i])
        print('\n')

        firstalignedseq = ""
        secondalignedseq = ""
        (I, J) = Index

        while tmatrix[I][J] != 0:
            if tmatrix[I][J] == 3:
                secondalignedseq += y[I - 1]
                firstalignedseq += x[J - 1]
                I -= 1
                J -= 1

            elif tmatrix[I][J] == 2:
                secondalignedseq += y[I - 1]
                firstalignedseq += '-'
                I -= 1

            else:
                secondalignedseq += '-'
                firstalignedseq += x[J - 1]
                J -= 1


        firstalignedseq = firstalignedseq[::-1]
        secondalignedseq = secondalignedseq[::-1]

        print(firstalignedseq)
        print(secondalignedseq)


matrix= open("blosum62.txt", 'r')
blosumfile = blosummatrix(matrix)

type_input=input('write type of sequence alignment:')
if type_input=='DNA' or type_input=='dna':
   firstsequence_input=input('write first sequence:')
   secondsequence_input=input('write second sequence:')
   LocalAlignmentforDNA(firstsequence_input, secondsequence_input)
elif type_input=='protein':
    firstsequence_input = input('write first sequence:')
    secondsequence_input = input('write second sequence:')
    LocalAlignmentforProtein(firstsequence_input, secondsequence_input,blosumfile)