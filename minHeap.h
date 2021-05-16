/*******************************************************************************
File minHeap.h is a header file for program KMC.c.
*******************************************************************************/

/* Function prototypes ********************************************************/
// returns the index of the parent node
int parent_node(int i) {
    return (i-1) / 2;
}

// return the index of the left child 
int left_child(int i) {
    return 2*i + 1;
}

// return the index of the right child 
int right_child(int i) {
    return 2*i + 2;
}
// swaps values of two integer variables
void swap_int(int *x, int *y) {
    int temp = *x;
    *x = *y;
    *y = temp;
}
// swaps values of two double variables
void swap_double(double *z, double *w) {
    double tmp = *z;
    *z = *w;
    *w = tmp;
}

// moves the item at position i of array a
// into its appropriate position
void min_heapify(int a[], double b[], int i, int n){
    // find left child node
    int left = left_child(i);

    // find right child node
    int right = right_child(i);

    // find the smallest among 3 nodes
    int smallest = i;

    // check if the left node is smaller than the current node
    if (left < n && b[left] < b[smallest]) {
        smallest = left;
    }

    // check if the right node is smaller than the current node
    if (right < n && b[right] < b[smallest]) {
        smallest = right;
    }

    // swap the smallest node with the current node 
    // and repeat this process until the current node is smaller than 
    // the right and the left node
    if (smallest != i) {
        swap_int(&a[i],&a[smallest]);
        swap_double(&b[i],&b[smallest]);

        min_heapify(a, b, smallest, n);
    }

}

void decrease_key(int a[], double b[], int i) {
    // find parent node
    int parent = parent_node(i);
    int largest = i;

    // check if the parent node is larger than the current node
    if(parent >= 0 && b[parent] > b[largest]) {
        largest = parent;
    }

    // swap the parent node with the current node if parent node is larger
    // and repeat this process until the current node is larger than 
    // the parent node
    if (largest != i) {
        swap_int(&a[i],&a[largest]);
        swap_double(&b[i],&b[largest]);
        
        decrease_key(a, b, largest);
    }
}

// converts an array into a binary min heap
void build_min_heap(int a[], double b[], int n) {
    for (int i = n/2; i >= 0; i--) {
        min_heapify(a, b, i, n);
    } 
}

// prints the heap
void print_heap(int a[], double b[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d\t%f\n", a[i],b[i]);
    }
    printf("\n");
}