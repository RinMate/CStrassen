/* strassen.c */
#include<stdio.h>
#include<stdlib.h>

#define FOR(i,a,b) for(int i=(a);i<(b);++i)
#define REP(i,n)   FOR(i,0,n)
#define REPR(i,n)  for(int i=(n);i>=(0);--i)

#define NEW(p,n) {p = malloc((n)*sizeof(p[0]));} 

// 行列の構造体
typedef struct {
	int n, m;   // n: 行数，m: 列数
	double* A;  // データ本体
}* dmatrix;

void dmatrix_strassen_init(dmatrix, dmatrix, dmatrix*, dmatrix*);
dmatrix dmatrix_strassen_recursive(dmatrix, dmatrix);

// 行列を生成する（全値を0.0で初期化していることに注意）
dmatrix dmatrix_new(int n, int m) {
	dmatrix M;
	NEW(M, 1);
	M->n = n;
	M->m = m;
	double* A;  // 本体
	NEW(A, m * n);
	REP(i, m * n) A[i] = 0.0;
	M->A = A;
	return M;
}

// 行列のメモリを解放
void dmatrix_free(dmatrix M){
	free(M->A);
	free(M);
}

// 行列の単位行列を返す
dmatrix dmatrix_identity(int n){
	dmatrix I = dmatrix_new(n, n);
	REP(i, n){
		I->A[i*n + i] = 1.0;
	}
	return I;
}

// N*N行列同士の積を返す（Strassenのアルゴリズム）
dmatrix dmatrix_strassen(dmatrix A, dmatrix B){
	if(A->n != A->m || A->n != B->n || B->n != B->m) return dmatrix_new(0, 0);

	dmatrix A_n, B_n;
	dmatrix_strassen_init(A, B, &A_n, &B_n);  // パディングして2のべき乗のサイズにする
	dmatrix C_n = dmatrix_strassen_recursive(A_n, B_n);
	dmatrix C = dmatrix_new(A->n, A->n);
	REP(i, A->n) REP(j, A->n) C->A[i*A->n + j] = C_n->A[i*A_n->n + j];

	dmatrix_free(C_n); dmatrix_free(B_n); dmatrix_free(A_n); 
	return C;
}

// 2のべき乗のサイズになるようパディング
void dmatrix_strassen_init(dmatrix A, dmatrix B, dmatrix *C, dmatrix *D){
	int n=1;
	while(n < A->n) n *= 2;
	dmatrix A_new = dmatrix_new(n, n);
	dmatrix B_new = dmatrix_new(n, n);
	REP(i, A->n){
		REP(j, A->n){
			A_new->A[i*n + j] = A->A[i*A->n + j];
			B_new->A[i*n + j] = B->A[i*A->n + j];
		}
	}
	*C = A_new;
	*D = B_new;
}

// 行列同士の和を返す
dmatrix dmatrix_add(dmatrix A, dmatrix B){
	if(A->n != B->n || A->m != B->m) return dmatrix_new(0, 0);

	int n = A->n;
	int m = A->m;
	dmatrix C = dmatrix_new(n, m);
	REP(i, n) REP(j, m) C->A[i*m + j] = A->A[i*m + j] + B->A[i*m + j];
	return C;
}

// 行列同士の差を返す
dmatrix dmatrix_sub(dmatrix A, dmatrix B){
	if(A->n != B->n || A->m != B->m) return dmatrix_new(0, 0);

	int n = A->n;
	int m = A->m;
	dmatrix C = dmatrix_new(n, m);
	REP(i, n) REP(j, m) C->A[i*m + j] = A->A[i*m + j] - B->A[i*m + j];
	return C;
}

// Strassenのアルゴリズム本体
dmatrix dmatrix_strassen_recursive(dmatrix A, dmatrix B){
	int N = A->n;
	if(N == 1){
		dmatrix C = dmatrix_new(1, 1);
		C->A[0] = A->A[0] * B->A[0];
		return C;
	}

	int n = N/2;  // Nは2のべき乗であることが保証されている

	dmatrix *As, *Bs;
	NEW(As, 4);
	NEW(Bs, 4);
	REP(i, 4){
		As[i] = dmatrix_new(n, n);
		Bs[i] = dmatrix_new(n, n);
	}

	REP(i, n) REP(j, n) As[0]->A[i*n + j] = A->A[i*N + j];
	REP(i, n) REP(j, n) As[1]->A[i*n + j] = A->A[i*N + (n+j)];
	REP(i, n) REP(j, n) As[2]->A[i*n + j] = A->A[(n+i)*N + j];
	REP(i, n) REP(j, n) As[3]->A[i*n + j] = A->A[(n+i)*N + (n+j)];
	REP(i, n) REP(j, n) Bs[0]->A[i*n + j] = B->A[i*N + j];
	REP(i, n) REP(j, n) Bs[1]->A[i*n + j] = B->A[i*N + (n+j)];
	REP(i, n) REP(j, n) Bs[2]->A[i*n + j] = B->A[(n+i)*N + j];
	REP(i, n) REP(j, n) Bs[3]->A[i*n + j] = B->A[(n+i)*N + (n+j)];

	dmatrix *tmps;
	NEW(tmps, 10);
	tmps[0] = dmatrix_add(As[0], As[3]);  // (A11+A22)
	tmps[1] = dmatrix_add(Bs[0], Bs[3]);  // (B11+B22)
	tmps[2] = dmatrix_add(As[2], As[3]);  // (A21+A22)
	tmps[3] = dmatrix_sub(Bs[1], Bs[3]);  // (B12-B22)
	tmps[4] = dmatrix_sub(Bs[2], Bs[0]);  // (B21-B11)
	tmps[5] = dmatrix_add(As[0], As[1]);  // (A11+A12)
	tmps[6] = dmatrix_sub(As[2], As[0]);  // (A21-A11)
	tmps[7] = dmatrix_add(Bs[0], Bs[1]);  // (B11+B12)
	tmps[8] = dmatrix_sub(As[1], As[3]);  // (A12-A22)
	tmps[9] = dmatrix_add(Bs[2], Bs[3]);  // (B21+B22)

	dmatrix *P;
	NEW(P, 7);
	P[0] = dmatrix_strassen_recursive(tmps[0], tmps[1]);  // P1 = (A11+A22)(B11+B22)
	P[1] = dmatrix_strassen_recursive(tmps[2], Bs[0]);  // P2 = (A21+A22)B11
	P[2] = dmatrix_strassen_recursive(As[0], tmps[3]);  // P3 = A11(B12-B22)
	P[3] = dmatrix_strassen_recursive(As[3], tmps[4]);  // P4 = A22(B21-B11)
	P[4] = dmatrix_strassen_recursive(tmps[5], Bs[3]);  // P5 = (A11+A12)B22
	P[5] = dmatrix_strassen_recursive(tmps[6], tmps[7]);  // P6 = (A21-A11)(B11+B12)
	P[6] = dmatrix_strassen_recursive(tmps[8], tmps[9]);  // P7 = (A12-A22)(B21+B22)

	dmatrix C = dmatrix_new(N, N);
	REP(i, n) REP(j, n) C->A[i*N + j] = P[0]->A[i*n + j] + P[3]->A[i*n + j] - P[4]->A[i*n + j] + P[6]->A[i*n + j];  // C11 = P1 + P4 - P5 + P7
	REP(i, n) REP(j, n) C->A[i*N + (n+j)] = P[2]->A[i*n + j] + P[4]->A[i*n + j];  // C12 = P3 + P5
	REP(i, n) REP(j, n) C->A[(n+i)*N + j] = P[1]->A[i*n + j] + P[3]->A[i*n + j];  // C21 = P2 + P4
	REP(i, n) REP(j, n) C->A[(n+i)*N + (n+j)] = P[0]->A[i*n + j] + P[2]->A[i*n + j] - P[1]->A[i*n + j] + P[5]->A[i*n + j];  // C22 = P1 + P3 - P2 + P6

	REPR(i, 3) { dmatrix_free(Bs[i]); dmatrix_free(As[i]); } 
	free(Bs); free(As); 
	REPR(i, 9) dmatrix_free(tmps[i]);
	free(tmps); 
	REPR(i, 6) dmatrix_free(P[i]);
	free(P); 

	return C;
}

// 行列のk乗を返す
dmatrix dmatrix_power(dmatrix A, int k){
	if(A->n != A->m) return dmatrix_new(0, 0);

	dmatrix R;
	if(k == 0) R = dmatrix_identity(A->n);
	else if(k == 1){
		dmatrix tmp = dmatrix_identity(A->n);
		R = dmatrix_strassen(A, tmp);
		dmatrix_free(tmp);
	}
	else {
		if(k%2 == 0){
			dmatrix tmp = dmatrix_power(A, k/2);
			R = dmatrix_strassen(tmp, tmp);
			dmatrix_free(tmp);
		}
		else{
			dmatrix tmp = dmatrix_power(A, (k-1)/2);
			R = dmatrix_strassen(A, dmatrix_strassen(tmp, tmp));
			dmatrix_free(tmp);
		}
	}
	return R;
}

// 行列を標準入力から読み取り
dmatrix dmatrix_read(){
	int n, m;
	scanf("%d %d", &n, &m);
	dmatrix M = dmatrix_new(n, m);

	double val;
	REP(i, n){
		REP(j, m){
			scanf("%lf", &val);
			M->A[i*m + j] = val;  
		}
	}

	return M;
}

// 行列を指定された形式で出力
void print_dmatrix(dmatrix M){
	int n = M->n;
	int m = M->m;
	printf("%d %d \n", n, m);

	double val;
	REP(i, n){
		REP(j, m){
			printf("%.10f ", M->A[i*m + j]);
		}
		printf("\n");
	}
}

int main(){
	dmatrix A, B;
	A = dmatrix_read();
	//print_dmatrix(A);

	int k;
	scanf("%d", &k);
	B = dmatrix_power(A, k);
	print_dmatrix(B);

	dmatrix_free(B);
	dmatrix_free(A);
	return 0;
}
