#ifndef MMA_FORMATTER_H_
#define MMA_FORMATTER_H_

/**
 * returns 'number' as a MMA compatible scientific notation string 
 * with 'precision' digits after the decimal point. The returned string
 * must be freed.
 * 
 * @brief converts double to a scientific notation string
 * @param number		the double to convert
 * @param precision	number of digits after decimal point in sci-not
 * @return 				the sci-not string e.g."-4.321*10^-03", which
 * 						should be freed
 */
char* getScientificNotation(double number, int precision);

/**
 * returns a string containing a MMA array of numbers in
 * scientific notation with `precision` digits after each decimal point.
 * The returned string must be freed.
 * 
 * @brief converts double array to a MMA array string of sci-not decimals
 * @param array			an array of doubles to convert to sci-not
 * @param length		number of doubles in array
 * @param precision	number of digits after decimal point in sci-not
 * @return 				a MMA array of sci-not e.g."{1.2*10^+3, ...}",
 * 						which should be freed
 */
char* convertDoubleArrToMMA(double* array, int length, int precision);

/**
 * returns a file handle to be passed to subsequent functions for
 * writing data to a MMA Association. File handle must be eventually
 * passed to closeAssocWrite.
 * 
 * @brief begins writing an association to file
 * @param filename		name of the file to save the MMA association to.
 * 						The assocation is read into MMA by Get[filename]
 * @return 				a file handle to pass to subsequent Association
 * 						Write functions, and which must be eventually
 * 						passed to closeAssocWrite
 */
FILE* openAssocWrite(char* filename);

/**
 * returns a file handle to be passed to subsequent functions for
 * writing data to a MMA Association, which appends elements to the
 * existing assoc in filename. File handle must be eventually
 * passed to closeAssocWrite or closeAssocAppend
 * 
 * @brief begins appending to an existing association in a file
 * @param filename		name of the file to containing the existing assoc
 * 						The assocation is read into MMA by Get[filename]
 * @return 				a file handle to pass to subsequent Association
 * 						Write functions, and which must be eventually
 * 						passed to closeAssocWrite or closeAssocAppend
 */
FILE* openAssocAppend(char* filename);

/**
 * @brief equivalent to closeAssocWrite
 */
void closeAssocAppend(FILE* file);

/**
 * @brief adds an int to the association
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param num		integer to add under key in association
 */
void writeIntToAssoc(FILE* file, char* keyname, int num);

/**
 * @brief adds a sci-not number to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param num			number to add in sci-not to the association
 * @param precision	number of digits after decimal point in sci-not
 */
void writeDoubleToAssoc(FILE* file, char* keyname, double num, int precision);

/**
 * @brief adds a string to the association
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param string	string to add under key in the association
 */
void writeStringToAssoc(FILE* file, char* keyname, char* string);

/**
 * @brief adds a MMA array of integers to the association
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param arr		array of ints to add under key
 * @param length	length of the array of ints
 */
void writeIntArrToAssoc(FILE* file, char* keyname, int* arr, int length);

/**
 * @brief adds a once nested array of integers, where the inner arrays have
 * an inconsistent length, to the association
 * @param file				file handle returned by openAssocWrite
 * @param keyname			key to add to the association
 * @param arr				2D array of ints to add under key
 * @param outerLength		number of inner arrays
 * @param innerLengths		length of each inner array
 * @param innerSpace		the allocated size of the inner arrays 
 */
void writeUnevenOnceNestedIntArrToAssoc(FILE* file, char* keyname, int (*arr)[], int outerLength, int* innerLengths, int innerSpace);
/**
 * @brief adds an array of long unsigned ints to the association 
 * 		  as a MMA array of integers
 * @param file		file handle returned by openAssocWrite
 * @param keyname	key to add to the association
 * @param arr		array of long unsigned ints to add under key
 * @param length	length of the array
 */
void writeUnsignedLongArrToAssoc(FILE* file, char* keyname, unsigned long *arr, int length);

/**
 * @brief adds a MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			array of doubles to convert to sci-not
 * @param length		length of the array
 * @param precision	number of digits after decimal in sci-not
 */
void writeDoubleArrToAssoc(FILE* file, char* keyname, double* arr, int length, int precision);

/**
 * @brief adds a once nested array of doubles, where the inner arrays have
 * an inconsistent length, to the association in scientific notation
 * @param file				file handle returned by openAssocWrite
 * @param keyname			key to add to the association
 * @param arr				2D array of doubles
 * @param outerLength		number of inner arrays in arr
 * @param innerLengths		length of each inner array
 * @param innerSpace		allocated space of each inner array
 * @param precision		number of digits after decimal in sci-not
 */
void writeUnevenOnceNestedDoubleArrToAssoc(
	FILE* file, char* keyname, double (*arr)[], int outerLength, int* innerLengths, int innerSpace, int precision);
	
/**
 * @brief adds a once-nested MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			nested array of doubles to convert to sci-not
 * @param outerLength	length of the outer list
 * @param innerLength	lenght of every inner list
 * @param precision	number of digits after decimal in sci-not
 */
void writeOnceNestedDoubleListToAssoc(
	FILE* file, char* keyname, double** arr, int outerLength, int innerLength, int precision
);

/**
 * @brief adds a numDimensions-nested MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			nested array of doubles to convert to sci-not
 * @param lengths		array of length of each dimension
 * @param precision	number of digits after decimal in sci-not
 */
void writeNestedDoubleListToAssoc(
	FILE* file, char* keyname, void* arr, int numDimensions, int* lengths, int precision
);

/**
 * @brief adds a numDimensions-nested MMA array of sci-not numbers to the association
 * @param file			file handle returned by openAssocWrite
 * @param keyname		key to add to the association
 * @param arr			nested array of doubles to convert to sci-not
 * @param lengths		array of length of each dimension
 * @param innerTrimLength	
 * 						how many of each inner-most array to keep (to trim out place-holder data)
 * @param precision	number of digits after decimal in sci-not
 */
void writeNestedDoubleArrToAssoc(
	FILE* file, char* keyname, void* arr, int numDimensions, int* lengths, int innerTrimLength, int precision
);

/**
 * @brief formats and finalises the association, so it is ready
 * 		   to read by MMA's Get[filename] function
 * @param file	file handle returned by openAssocWrite
 */
void closeAssocWrite(FILE* file);

/**
 * @brief opens an existing MMA association (which must not be empty)
 *        so that new keys may be appended to it. Note writing new
 *        data to an existing key element will cause duplication in-file,
 *        and the new data will be ignored by Mathematica when parsed
 * @return 				a file handle to pass to subsequent Association
 * 						Write functions, and which must be eventually
 * 						passed to closeAssocWrite or closeAssocAppend
 */
FILE* openAssocAppend(char* filename);

/**
 * @brief formats and finalises the association, so it is ready
 * 		   to read by MMA's Get[filename] function
 * @param file	file handle returned by openAssocAppend
 */
void closeAssocAppend(FILE* file);
	
#endif // MMA_FORMATTER_H_

