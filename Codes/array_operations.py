#Array operations

import math
#finding position of a value in array
def find_index(X,val):
    n=len(X)
    for i in range(n):
        if round(X[i],1)==round(val,1):
            return(i)

        #
    #


# Returns index of x in arr if present, else -1
def binarySearch(arr, l, r, x):
    # Check base case
    if r >= l:

        mid = l + (r - l) // 2

        # If element is present at the middle itself
        if round(arr[mid],1) == round(x,1):
            return mid

            # If element is smaller than mid, then it
        # can only be present in left subarray
        elif round(arr[mid],1) > round(x,1):
            return binarySearch(arr, l, mid - 1, x)

            # Else the element can only be present
        # in right subarray
        else:
            return binarySearch(arr, mid + 1, r, x)

