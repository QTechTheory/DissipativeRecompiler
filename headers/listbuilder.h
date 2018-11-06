#ifndef LIST_BUILDER_H_
#define LIST_BUILDER_H_


double** allocDoublyNestedDoubleList(int outerLen, int innerLen);

double*** allocTriplyNestedDoubleList(int outerLen, int middleLen, int innerLen);

void freeDoublyNestedDoubleList(double** list, int outerLen);

void freeTriplyNestedDoubleList(double*** list, int outerLen, int middleLen);

#endif // LIST_BUILDER_H_