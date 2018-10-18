#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "bmp24_lib.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>

bmp24_file::bmp24_file(WORD width, WORD height)
{
	common_init(width, height);
}

bmp24_file::bmp24_file(WORD width, WORD height, std::string filename)
{
	common_init(width, height);
	name = filename;
}

bmp24_file::~bmp24_file()
{
	for (unsigned int i = 0; i < size_height; i++)
		delete[] color_array[i];
	delete[] color_array;
}

void bmp24_file::common_init(WORD width, WORD height)
{
	size_width = width; // высота
	size_height = height; // ширина
	color_array = new RGBTRIPLE *[height];
	// выделение памяти под все пиксели?
	for (unsigned int i = 0; i < height; i++)
	{
		color_array[i] = new RGBTRIPLE[width];
		memset(color_array[i], 0, sizeof(RGBTRIPLE) * width);
	}

	memset(&file_header, 0, sizeof(BITMAPFILEHEADER));
	// Тип файла (BM)
	file_header.bfType = 0x4d42;
	// Смещение от начала битов 
	file_header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	// Размер файла
	file_header.bfSize = file_header.bfOffBits + sizeof(RGBTRIPLE) * size_width * size_height + size_height * ((sizeof(RGBTRIPLE) * size_width) % 4);

	memset(&info_header, 0, sizeof(BITMAPINFOHEADER)); // так положено, лол
	info_header.biSize = sizeof(BITMAPINFOHEADER);
	info_header.biWidth = size_width;
	info_header.biHeight = size_height;
	info_header.biPlanes = 1;
	info_header.biBitCount = 24; // 24 бит на пиксель
	info_header.biCompression = 0;

	// Выделяем память под пустой цвет
	memset(&null_color, 0, sizeof(RGBTRIPLE));
}


void bmp24_file::set_vertical_mirror(bool value)
{
	//  если value = true и начало координат - левый нижний угол или наоборот
	if ((value && info_header.biHeight > 0) || (!value && info_header.biHeight < 0))
		info_header.biHeight *= -1; // меняем начало координат на противоположное
}

void bmp24_file::set_pixel(WORD x, WORD y, RGBTRIPLE c)
{	// если выходит за рамки точечного рисунка, не озвращаем значения
	if (x >= size_width || y >= size_height)
		return;
	// устанавливаем цвет пикселя
	color_array[y][x] = c; // почму наоборот?
}

void bmp24_file::set_pixel(WORD x, WORD y, BYTE r, BYTE g, BYTE b)
{
	RGBTRIPLE c;
	// устанавливам компоненты цветов
	c.rgbtBlue = b;
	c.rgbtGreen = g;
	c.rgbtRed = r;
	set_pixel(x, y, c);
}

RGBTRIPLE bmp24_file::get_pixel(WORD x, WORD y)
{
	// если выходим за рамки пиксельного рисунка, возвращаем пустой цвет
	if (x >= size_width || y >= size_height)
		return null_color;
	return color_array[y][x]; // возвращаем цветовое значение пикселя
}

void bmp24_file::get_pixel(WORD x, WORD y, RGBTRIPLE & c)
{
	c = get_pixel(x, y);
}

void bmp24_file::get_pixel(WORD x, WORD y, BYTE & r, BYTE & g, BYTE & b)
{
	RGBTRIPLE c = get_pixel(x, y);
	r = c.rgbtRed;
	g = c.rgbtGreen;
	b = c.rgbtBlue;
}

void bmp24_file::write()
{
	FILE * f = fopen(name.c_str(), "wb");
	if (!f)
	{
		if (name.length() == 0)
			std::cerr << "Error: need to set name of file!" << std::endl;
		else
			std::cerr << "Error: unable to open file \"" << name << "\"!" << std::endl;
		return;
	}
	/*	size_t fwrite(const void *buf, size_t size, size_t count, FILE *stream)
		Функция fwrite() записывает count объектов — каждый объект по size символов в длину — в поток,
		указанный stream, из символьного массива, указанного buf. Указатель позиции в файле продвигается 
		вперед на количество записанных символов.
		Функция fwrite() возвращает количество действительно записанных объектов, которое в случае успеха 
		равно затребованному количеству. Если количество записанных объектов меньше, чем это указано при
		вызове, то произошла ошибка.
	*/
	fwrite(&file_header, sizeof(BITMAPFILEHEADER), 1, f);
	fwrite(&info_header, sizeof(BITMAPINFOHEADER), 1, f);
	for (unsigned int i = 0; i < size_height; i++)
	{
		fwrite(color_array[i], sizeof(RGBTRIPLE), size_width, f);
		fwrite(&null_color, (sizeof(RGBTRIPLE) * size_width) % 4, 1, f);
	}
	fflush(f);
	fclose(f);
}

void bmp24_file::write(std::string filename)
{
	name = filename;
	write();
}
