# Lap-trinh-song-song

Phương trình nhiệt 2 chiều

Phương pháp Gauss Seidel

1. THAM SỐ ĐẦU VÀO

#define m 20

#define n 20

#define epsilon    0.001 

#define tolerance  0.001


2. KHỞI TẠO ĐIỀU KIỆN BAN ĐẦU


void KhoiTao(float *C) {
	
	int i, j;
	
	for (int i = 0; i < m; i++)
	
		for (int j = 0; j < n; j++)
		
			if (i >= (m / 2 - 5) && i < (m / 2 + 5) 
			
			 && j >= (n / 2 - 5) && j < (m / 2 + 5))
			 
				*(C + n * i + j) = 80;
				
			else
			
				*(C + n * i + j) = 25;				
}
