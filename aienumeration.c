#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_permutation.h>

#define LEVEL 3
#define BLOCKS 3
#define N 9

// Function prototypes
gsl_matrix* CreateBlockMatrix(int rows, int cols, int* blockshape, gsl_matrix** blocklist);
int is_zeroone(gsl_matrix* matrix);
gsl_matrix* switchingblock(int level, int blocknumber, const char* sort);
int* getvertices(int n, int* S, int len_S);
void symmetrymatrices(gsl_matrix*** matrices, int* count, int n, int blocks, int level, const char* sort, int recursionlevel);
void Bruteswitching(gsl_matrix** gen, int n, int level, int blocks, const char* sort);

// Function to create a block matrix
gsl_matrix* CreateBlockMatrix(int rows, int cols, int* blockshape, gsl_matrix** blocklist) {
    gsl_matrix* result = gsl_matrix_alloc(rows * LEVEL, cols * LEVEL);
    gsl_matrix_set_zero(result);

    int block_index = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            gsl_matrix* block = blocklist[blockshape[block_index++]];
            for (int bi = 0; bi < block->size1; bi++) {
                for (int bj = 0; bj < block->size2; bj++) {
                    gsl_matrix_set(result, i * LEVEL + bi, j * LEVEL + bj, gsl_matrix_get(block, bi, bj));
                }
            }
        }
    }

    return result;
}

// Function to check if a matrix is a zero-one matrix
int is_zeroone(gsl_matrix* matrix) {
    for (size_t i = 0; i < matrix->size1; i++) {
        for (size_t j = 0; j < matrix->size2; j++) {
            double value = gsl_matrix_get(matrix, i, j);
            if (value != 0.0 && value != 1.0) {
                return 0;  // Not a zero-one matrix
            }
        }
    }
    return 1;
}

// Function to create a switching block matrix
gsl_matrix* switchingblock(int level, int blocknumber, const char* sort) {
    gsl_matrix* J = gsl_matrix_alloc(level, level);
    gsl_matrix* O = gsl_matrix_alloc(level, level);
    gsl_matrix* Y = gsl_matrix_alloc(level, level);

    // Initialize matrices
    gsl_matrix_set_all(J, 1.0);
    gsl_matrix_set_zero(O);
    gsl_matrix_set_identity(Y);
    gsl_matrix_scale(Y, level);
    gsl_matrix_sub(Y, J);

    int baserow[3] = {0, 2, 1};  // Assumes blocknumber >= 2
    int shapematrix[blocknumber * blocknumber];
    int index = 0;

    for (int i = 0; i < blocknumber; i++) {
        for (int j = 0; j < blocknumber; j++) {
            shapematrix[index++] = baserow[j % 3];
        }
        int temp = baserow[2];
        baserow[2] = baserow[1];
        baserow[1] = baserow[0];
        baserow[0] = temp;
    }

    gsl_matrix* result;
    if (strcmp(sort, "wqh") == 0) {
        gsl_matrix* blocks[] = {Y, O, J};
        result = CreateBlockMatrix(blocknumber, blocknumber, shapematrix, blocks);
        gsl_matrix_scale(result, 1.0 / level);
    } else if (strcmp(sort, "gm") == 0) {
        gsl_matrix* blocks[] = {Y, O, J};
        result = CreateBlockMatrix(blocknumber, blocknumber, shapematrix, blocks);
        gsl_matrix_scale(result, -1.0 / level);
    } else {
        result = NULL;
    }

    gsl_matrix_free(J);
    gsl_matrix_free(O);
    gsl_matrix_free(Y);

    return result;
}

// Function to get vertices with a permutation
int* getvertices(int n, int* S, int len_S) {
    int* vertices = (int*)malloc(n * sizeof(int));
    int k = 0;

    for (int i = 0; i < len_S; i++) {
        vertices[k++] = S[i];
    }

    for (int i = 1; i <= n; i++) {
        int found = 0;
        for (int j = 0; j < len_S; j++) {
            if (S[j] == i) {
                found = 1;
                break;
            }
        }
        if (!found) {
            vertices[k++] = i;
        }
    }

    return vertices;
}

// Function to generate symmetry matrices recursively
void symmetrymatrices(gsl_matrix*** matrices, int* count, int n, int blocks, int level, const char* sort, int recursionlevel) {
    if (recursionlevel == 1) {
        gsl_matrix* identity = gsl_matrix_alloc(n, n);
        gsl_matrix_set_identity(identity);
        matrices[*count] = identity;
        (*count)++;
        return;
    }

    gsl_combination* c = gsl_combination_alloc(n, level);
    gsl_combination_init_first(c);

    do {
        int S[level];
        for (int i = 0; i < level; i++) {
            S[i] = gsl_combination_get(c, i) + 1;
        }

        int* perm = getvertices(n, S, level);
        gsl_permutation* p = gsl_permutation_alloc(n);
        gsl_permutation_init(p);

        gsl_matrix* Xmatrix = gsl_matrix_alloc(n, n);
        for (int i = 0; i < n; i++) {
            gsl_permutation_set(p, i, perm[i] - 1);
        }
        gsl_permute_matrix(p, Xmatrix);

        int current_count = *count;
        symmetrymatrices(matrices, count, n, blocks, level, sort, recursionlevel - 1);

        for (int i = current_count; i < *count; i++) {
            gsl_matrix* P = matrices[i];
            gsl_matrix* result = gsl_matrix_alloc(n, n);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Xmatrix, P, 0.0, result);
            matrices[*count] = result;
            (*count)++;
        }

        free(perm);
        gsl_permutation_free(p);
        gsl_matrix_free(Xmatrix);
    } while (gsl_combination_next(c) == GSL_SUCCESS);

    gsl_combination_free(c);
}

// Brute-force switching function
void Bruteswitching(gsl_matrix** gen, int n, int level, int blocks, const char* sort) {
    gsl_matrix* Q = switchingblock(level, blocks, sort);
    gsl_matrix* identity = gsl_matrix_alloc(n - blocks * level, n - blocks * level);
    gsl_matrix_set_identity(identity);

    gsl_matrix* Q_extended = gsl_matrix_alloc(n, n);
    gsl_matrix_view Q_view = gsl_matrix_submatrix(Q_extended, 0, 0, Q->size1, Q->size2);
    gsl_matrix_memcpy(&Q_view.matrix, Q);

    gsl_matrix* result = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(result);
    gsl_matrix_add(result, Q_extended);
    gsl_matrix_free(Q_extended);

    int count_total = 0;
    int count_switching_sets = 0;

    int num_permutations = 0;
    gsl_matrix** perm_matrices = malloc(1000 * sizeof(gsl_matrix*));
    symmetrymatrices(perm_matrices, &num_permutations, n, blocks, level, sort, blocks);

    for (int g = 0; g < 100; g++) {  // Placeholder loop, assuming gen has 100 graphs
        gsl_matrix* A = gen[g];  // Assume gen[g] gives the adjacency matrix of the g-th graph
        if (A == NULL) continue;

        int count_cospectral = 0;
        count_total++;

        gsl_combination* combinations = gsl_combination_alloc(n, blocks * level);
        gsl_combination_init_first(combinations);

        do {
            int S[blocks * level];
            for (int i = 0; i < blocks * level; i++) {
                S[i] = gsl_combination_get(combinations, i) + 1;
            }

            int* perm = getvertices(n, S, blocks * level);
            gsl_permutation* p0 = gsl_permutation_alloc(n);
            for (int i = 0; i < n; i++) {
                gsl_permutation_set(p0, i, perm[i] - 1);
            }

            gsl_matrix* P0 = gsl_matrix_alloc(n, n);
            gsl_permute_matrix(p0, P0);

            for (int i = 0; i < num_permutations; i++) {
                gsl_matrix* P1 = perm_matrices[i];
                gsl_matrix* P = gsl_matrix_alloc(n, n);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P0, P1, 0.0, P);

                gsl_matrix* B = gsl_matrix_alloc(n, n);
                gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, result, A, 0.0, B);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, B, P, 0.0, B);
                gsl_matrix_transpose(B);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, B, result, 0.0, B);

                if (is_zeroone(B)) {
                    count_switching_sets++;
                    // Check if B is isomorphic to A (this requires additional implementation)
                    // if (!Graph_is_isomorphic(B, A)) {
                    //     count_cospectral++;
                    // }
                }

                gsl_matrix_free(B);
                gsl_matrix_free(P);
            }

            free(perm);
            gsl_permutation_free(p0);
            gsl_matrix_free(P0);
        } while (gsl_combination_next(combinations) == GSL_SUCCESS);

        gsl_combination_free(combinations);

        // if (count_cospectral > 0) {
        //     Add G to switchgraphs
        //     Add count_cospectral to counts
        // }
    }

    gsl_matrix_free(identity);
    gsl_matrix_free(Q);
    for (int i = 0; i < num_permutations; i++) {
        gsl_matrix_free(perm_matrices[i]);
    }
    free(perm_matrices);
}