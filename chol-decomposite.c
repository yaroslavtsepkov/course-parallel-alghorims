#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
double **alloc_matrix(int row, int column){
    double **matrix;
    matrix=malloc(row*sizeof(double*));
    matrix[0] = (double *)malloc(row*column*sizeof(double));
    for(int i = 1; i < row; i++)
        matrix[i] = matrix[0]+i*column;
    return matrix;
}
double **eye(int row, int column){
    double **eye_matrix;
    for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            if(i==j){
                eye_matrix[i][j] = 1;
            }else{
                eye_matrix[i][j] = 0;
            }
        }
    }
    return eye_matrix;
}
// void fill(int row, int column, double matrix[row][column]){
//     for(int i=0; i<row; i++){
//         for(int j=0; j<column; j++)
//             matrix[i][j]=rand()%10;
//     }
// }
void show(int row, int column, double matrix[row][column]){
    for(int i=0; i<row; i++){
        for(int j=0; j<column; j++)
            printf("%.2lf ", matrix[i][j]);
        printf("\n");
    }
}
// void chol_bana_decomp(int row, int column, double matrix[row][column]){
//     for (int i = 0; i < row; i++)
//     {
//         for (int j = 0; j < i; j++)
//         {
//             double sum = 0;
//             for (int k = 0; k < j; k++){
//                 sum+= matrix[i][k]*matrix[j][k];
//             if (i == j)
//                 matrix[i][j] = sqrt(matrix[i][j]-sum);
//             else
//                 matrix[i][j] = (1.0 / matrix[j][j]*(matrix[i][j]-matrix[i][k]*matrix[j][k]));
//             }
//         }
//     }
    
// }
// void chol_crout_decomp(int row, int column, double matrix[row][column]){
//     for (int j = 0; j < row; j++) {
//         float sum = 0;
//         for (int k = 0; k < j; k++) {
//             sum += matrix[j][k] * matrix[j][k];
//         }
//         matrix[j][j] = sqrt(matrix[j][j] - sum);

//         for (int i = j + 1; i < row; i++) {
//             sum = 0;
//             for (int k = 0; k < j; k++) {
//                 sum += matrix[i][k] * matrix[j][k];
//             }
//             matrix[i][j] = (1.0 / matrix[j][j] * (matrix[i][j] - sum));
//         }
//     }
// }
// void preprocessing(int row, int column, double matrix[row][column]){
//     double **temp;
//     double **eye_matrix;
//     temp = alloc_matrix(row, column);
//     eye(row, column, (double(*)[])(*eye_matrix));
//     show(row, column, (double(*)[])(*eye_matrix));
//     for (int i = 0; i < row; i++)
//     {
//         for (int j = 0; j < column; j++)
//         {
//             temp[i][j] = matrix[j][i];
//         }
        
//     }
//      for (int i = 0; i < row; i++)
//     {
//         for (int j = 0; j < column; j++)
//         {
//             matrix[i][j] = 0.5*(matrix[i][j] + temp[i][j]);
//         }
        
//     }
//     free(*eye_matrix);
//     free(eye_matrix);
//     free(*temp);
//     free(temp);
// }
int main(int argc, char const *argv[])
{
    int row=3, col=3;
    double **eye_matrix;
    printf("eye matrix\n");
    eye_matrix = eye(row, col);
    show(row,col,(double(*)[])(*eye_matrix));
    // int row,column;
    // double **A;
    // printf("enter rows: ");
    // scanf("%d", &row);
    // printf("enter columns: ");
    // scanf("%d", &column);
    // A = alloc_matrix(row, column);
    // fill(row, column, (double (*)[])(*A));
    // printf("%s \n", "source matrix");
    // show(row, column, (double (*)[])(*A));
    // printf("%s \n", "after preprocessing");
    // preprocessing(row, column, (double (*)[])(*A));
    // show(row, column, (double (*)[])(*A));
    // printf("%s \n", "cholesky decomposit");
    // chol_crout_decomp(row, column,(double (*)[])(*A));
    // show(row, column, (double (*)[])(*A));
    // free(*A);
    // free(A);
    return 0;
}
