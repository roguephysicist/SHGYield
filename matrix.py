import types
import numpy

from sys import argv

script, xyzfile = argv

# This function is to load an ascii format matrix (float numbers separated by whitespace characters and newlines) into a numpy matrix object. f is a file object or a file path. 
def load_matrix_from_file(f):
    if type(f) == types.StringType:
        fo = open(f, 'r')
        matrix = load_matrix_from_file(fo)
        fo.close()
        return matrix
    elif type(f) == types.FileType:
        file_content = f.read().strip()
        file_content = file_content.replace('\r\n', ';')
        file_content = file_content.replace('\n', ';')
        file_content = file_content.replace('\r', ';')
    return numpy.matrix(file_content)

    raise TypeError('f must be a file object or a file name.')

a = load_matrix_from_file(xyzfile)
print "The matrix to be converted is \n%s" % a
b = a*(0.529177249)
print "\nOur new matrix is now \n%s" % b
