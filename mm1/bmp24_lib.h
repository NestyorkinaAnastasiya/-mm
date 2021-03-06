#pragma once
#ifndef BMP24_LIB_H_INCLUDED
#define BMP24_LIB_H_INCLUDED

#include <string>

#ifdef _WIN32

#include <windows.h>

#else

#include <stdint.h>
typedef uint8_t  BYTE;
typedef uint16_t WORD;
typedef uint32_t DWORD;
typedef int32_t  LONG;

#pragma pack(push, 1)

typedef struct tagRGBTRIPLE
{
	BYTE rgbtBlue;
	BYTE rgbtGreen;
	BYTE rgbtRed;
} RGBTRIPLE;

typedef struct tagBITMAPFILEHEADER
{
	WORD bfType;
	DWORD bfSize;
	WORD bfReserved1;
	WORD bfReserved2;
	DWORD bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER
{
	DWORD biSize;
	LONG biWidth;
	LONG biHeight;
	WORD biPlanes;
	WORD biBitCount;
	DWORD biCompression;
	DWORD biSizeImage;
	LONG biXPelsPerMeter;
	LONG biYPelsPerMeter;
	DWORD biClrUsed;
	DWORD biClrImportant;
} BITMAPINFOHEADER;

#pragma pack(pop)

#endif
class bmp24_file
{
public:
	bmp24_file(WORD width, WORD height);
	bmp24_file(WORD width, WORD height, std::string filename);
	~bmp24_file();
	void set_vertical_mirror(bool value);
	void write();
	void write(std::string filename);
	void set_pixel(WORD x, WORD y, RGBTRIPLE c);
	void set_pixel(WORD x, WORD y, BYTE r, BYTE g, BYTE b);
	RGBTRIPLE get_pixel(WORD x, WORD y);
	void get_pixel(WORD x, WORD y, RGBTRIPLE & c);
	void get_pixel(WORD x, WORD y, BYTE & r, BYTE & g, BYTE & b);

protected:
	BITMAPFILEHEADER file_header;
	BITMAPINFOHEADER info_header;
	WORD size_width;
	WORD size_height;
	std::string name;
	RGBTRIPLE ** color_array;
	RGBTRIPLE null_color;
	void common_init(WORD width, WORD height);
};

#endif // BMP24_LIB_H_INCLUDED

/* RGBTRIPLE - ���������, ����������� ���� { rgbtBlue, rgbtGreen, rgbtRed }
 BITMAPFILEHEADER - ���������, ���������� ���. � ����, ������� � ������ �����,
				������� �������� ���������-����������� ��������� ������ DIB
				����� ���������: {	bfType - ��� ����� (must be BM :c)
									bfSize - ���������� ������ ����� ��������� ������� � ������
									bfReserved1, bfReserved2 - ���������������, ������ ���� 0
									bfOffBits - ���������� �������� �� ������ ����� ��������� �������
										� ��������� BITMAPFILEHEADER, � ������
								 }
 BITMAPINFOHEADER - ���������, ���������� ���. � �������� � �������� �������
					  ���������-������������ ��������� ������� DIB
				����� ���������: {	bfSize	 -	���������� ����� ������ ����������� ��� ���������
									bfWidth	 -	���������� ������ ��������� ������� � ��������.
												���� ��������� biWidth ������������� ������ ������������ ����������� �� ����� ����������� ������� JPEG ��� PNG , ��������������.
									biHeight -	������������� ������ ��������� �������, � ��������.���� ������������� ��������, �������� ������� - ���������� ���������-����������� ��������� ������(DIB),
												� ��� ������ ��������� - ����� ������ ����. ���� ������������� ��������, �������� ������� - ���������� ���������-����������� ��������� ������(DIB),
												� ��� ������ ��������� - ����� ������� ����. ���� biHeight ����� ������������� ��������, ������� ��������� �� ���������� ������(DIB), �������� biCompression
												������ ���� ��� BI_RGB ��� BI_BITFIELDS. ���������� ��������� �������(DIB) �� ����� ���������.
									biPlanes -  ������������� ����� ���������� �������� ����������.
									biBitCount - ������������� ����� ����� �� �������. ������������� ����� �����, ������� ���������� ������ ������� � ������������ ����� ������ � �������� �������.
													���� ���� ��������� ������ ���� ����� �� ����������������� ��������.*/
/*													0	Windows 98 / Me, Windows 2000 / XP: ����� ����� �� ������� �������� ��� ��������������� �������� ����������� PNG ��� JPEG.
													1	�������� ������� �������� �����������, � ���� bmiColors  ��������� BITMAPINFO �������� ��� ������. ������ ��� � ��������� ������� ������������ �������.
														���� ��� �������, ������� ������������ ������ ������ ������ � ������� bmiColors; ���� ��� ����������, ������� ����� ���� ������ ������ � �������.
													4	�������� ������� ����� �������� 16 ������, � ����  bmiColors ��������� BITMAPINFO �������� �� 16 �������.������ ������� � �������� ������� �����������
														4 - ������ �������� � ������� ������.��������, ���� ������ ���� � �������� ������� - 0x1F, �� ������������ ��� �������. ������ ������� �������� ���� ��
														������ ������ �������, � ������ ������� �������� ���� � ������������ ������ �������.
													8	�������� ������� ����� �������� 256 ������, � ����  bmiColors ��������� BITMAPINFO �������� �� 256 �������.� ���� ������, ������ ���� � ������� ������������
														��������� �������.
													16	�������� ������� ����� �������� 216 ������.���� ���� biCompression  ��������� BITMAPINFOHEADER ����� BI_RGB, �� ���� bmiColors  ��������� BITMAPINFO �����
														�������� �����(NULL).������ �����(WORD) � ������� ��������� ������� ������������ ��������� �������.������������� ������������� ��������, �������� � ������
														����� ������������ � ����� ������ ��� ������� ���������� �����.�������� ��� ������ ��������� � ����� ������� ���� �����, ��������� ���� ����� ��� �������
														�������� � �������� �����.������� �������� ��� �� ������������.������� ������ bmiColors ������������ ��� ����, ����� �������������� �����, ������������ ��
														����������� �������� �������, � ������ ��������� ����� �������, �������� ������ biClrUsed  ��������� BITMAPINFOHEADER.
													���� ���� biCompression  ��������� BITMAPINFOHEADER - BI_BITFIELDS, ���� ��������� bmiColors ��������  ��� ������� �����(DWORD) ����� �����, �������,
													��������������, ���������� �������, ������� � ����� ����������  ������� �������.������ �����(WORD) � ������� ��������� ������� ������������ ��������� �������.
													24	�������� ������� ����� �������� 2 24 ������, � ���� bmiColors ��������� BITMAPINFO ����� �������� �����(NULL).������ 3 - �������� ������� � ������� ���������
														������� ������������ ������������� ������������� ������, �������� � �������� �����, ��������������, ��� �������.������� ������ bmiColors  ������������ ��� ����,
														����� �������������� �����, ������������ �� ����������� �������� ������� � ������ ��������� ����� �������, ������������ ������ biClrUsed  ��������� BITMAPINFOHEADER.
													32	�������� ������� ����� �������� 232 ������.���� ���� biCompression  ��������� BITMAPINFOHEADER ����� BI_RGB, ���� bmiColors  ��������� BITMAPINFO ����� ��������
														�����(NULL).������ ������� �����(DWORD) � ������� ��������� ������� ������������ ������������� ������������� ������, �������� � �������� �����, ��������������, ���
														�������.������� ���� � ������ ������� �����(DWORD) �� ������������.������� ������ bmiColors  ������������ ��� ����, ����� �������������� �����, ������������ ��
														����������� �������� ������� � ������ ��������� ����� �������, ������������ ������ biClrUsed  ��������� BITMAPINFOHEADER. */
/*									biCompression - ���������� ��� ������ ��� ������� ������� ����� ����� ��������� �������(������ ������ ���� ��������� - ����������� ��������� ������(DIB) �� ����� ���������).
													���� ���� ��������� ����� ���� ����� �� ����������������� ��������.*/
/*													BI_RGB			����������� ������.
													BI_RLE8			������ � ������������ ����� �����(RLE) ��� �������� �������� � 8 ������ �� �������(bpp).���� ������ ������ - 2 - �������� ������, ��������� �� �����
																	����� ���������, ��������������� ������, ����������  ������ �����.��������� ���������� ��.� ������ ������ ��������� �������.
													BI_RLE4	RLE		������ ��� �������� �������� � 4 ������ �� �������(bpp).���� ������ ������ - 2 - �������� ������, ��������� �� ����� ����� ���������, ��������������
																	����� ��������� ����� ������ � �����(DWORD).��������� ���������� ��.� ������ ������ ��������� �������.
													BI_BITFIELDS	����������, ��� �������� ������� �� ���� �, ��� ������� ������ ������� �� ���� ������� ����(DWORD) ����� �����, ������� ������ �������, ������� � ����� ���������� �����, ��������������, ��� ������� �������.��� ��������� �����, ����� ������������ �������� ������� � 16 - � 32 - ������ �� �������(bpp).
													BI_JPEG	Windows 98 / Me, Windows 2000 / XP: ��������� �� ��, ��� ����������� ����� ������ JPEG.
													BI_PNG	Windows 98 / Me, Windows 2000 / XP : ��������� �� ��, ��� ����������� ����� ������� PNG.*/
/*									biSizeImage - ������������� ������ �����������, � ������.�� ����� ���� ���������� � ���� ��� BI_RGB �������� ��������.
									biXPelsPerMeter - ������������� ����������� ����������� �� ��������� ��� ������������ ����������  ��������� �������, � �������� �� ����. ���������� ����� ������������ ���
													  ��������, ����� �������� ������� ����� �� ������ ��������, ������� ����� ����� ������������� ��������������� �������� ����������.
									biYPelsPerMeter - ������������� ����������� ����������� �� ��������� ��� ������������ ����������  ��������� �������, � �������� �� ����.
									biClrUsed - ������������� �����  �������� ����� � ������� ������, ������� ���������� ������������ �������� ��������.���� ��� �������� ��������� ����, �������� �������
												���������� ������������ ����� ������, ��������������� �������� ����� ��������� biBitCount ���  ������ ������, ��������� ������ biCompression.
												���� biClrUsed �� ����, � ���� ��������� biBitCount - ������ ��� 16, ���� ��������� biClrUsed  ������������� ����������� ����� ������ ��������� �����������
												������ ��� �������� ����������.���� biBitCount  ����� 16 ��� ������, ���� ��������� biClrUsed  ������������� ������ ������������ �������� ������ ��� �����������
												��������� �������� ������.���� biBitCount ��������� 16 ��� 32, ����������� �������� ������� ������� ��������������� �� ����� �������� �������(DWORD) �����.
												����� ������ ��������� ������� ������� ��������������� �� ���������� BITMAPINFO, ��� - ������ �������� �������.�� ������ ������� ������ ��������� ��������� ����������.
												������ �������� ������� �������, ����� ���� ��������� biClrUsed ��� ��� ����� ��� ����������� �������� ������� ������.
									biClrImportant - ������������� ����� �������� �����, ������� ��������� ����� �������� �� ������ �������� �������.���� ��� �������� ��������� ����, ��������� ��� ����� .*/
