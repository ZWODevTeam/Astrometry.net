#include <iostream>
#include <numeric>
#include <chrono>
#include <opencv2/opencv.hpp>
#include "fitsio.h"

using std::cout;
using std::endl;

class Edge
{
public:
	virtual ~Edge() = default;
	float x{ 0 };
	float y{ 0 };
	int val{ 0 };
	int scanned{ 0 };
	float width{ 0 };
	float HFR{ -1 };
	float sum{ 0 };
};

inline unsigned	short round_to_WORD(double x) {
	if (x <= 0.0)
		return (unsigned short)0;
	if (x > USHRT_MAX)
		return USHRT_MAX;
	return (unsigned short)(x + 0.5);
}

void convert_data(int bitpix, const void *from, unsigned short *to, unsigned int nbdata, bool values_above_1) {
	int i;
	unsigned char *data8;
	int16_t *data16;
	double *pixels_double;
	double norm = 1.0;
	long *sdata32;	// TO BE TESTED on 32-bit arch, seems to be a cfitsio bug
	unsigned long *data32;

	switch (bitpix) {
	case BYTE_IMG:
		data8 = (unsigned char *)from;
		for (i = 0; i < nbdata; i++)
			to[i] = (unsigned short)data8[i];
		break;
	case USHORT_IMG:	// siril 0.9 native
						// nothing to do
		break;
	case SHORT_IMG:
		// add 2^15 to the read data to obtain unsigned
		data16 = (int16_t *)from;
		for (i = 0; i < nbdata; i++) {
			int sum = 32768 + (int)data16[i];
			to[i] = (unsigned short)sum;
		}
		break;
	case ULONG_IMG:		// 32-bit unsigned integer pixels
		data32 = (unsigned long *)from;
		for (i = 0; i < nbdata; i++)
			to[i] = (unsigned short)(data32[i] >> 16);
		break;
	case LONG_IMG:		// 32-bit signed integer pixels
		sdata32 = (long *)from;
		for (i = 0; i < nbdata; i++)
			to[i] = (unsigned short)((sdata32[i] >> 16) + 32768);
		break;
	case DOUBLE_IMG:	// 64-bit floating point pixels
	case FLOAT_IMG:		// 32-bit floating point pixels
		pixels_double = (double *)from;
		/* various data values can be found in a float or
		* double image. Sometimes it's normalized between 0
		* and 1, but sometimes, DATAMIN and DATAMAX give the range.
		*/
		if (!values_above_1) norm = USHRT_MAX;
		for (i = 0; i < nbdata; i++) {
			to[i] = round_to_WORD(norm * pixels_double[i]);
		}
		break;
	case LONGLONG_IMG:	// 64-bit integer pixels
	default:
		cout<<"Unknown FITS data format in internal conversion\n";
	}
}

int readfits(const char* filename, unsigned short** img, int& width, int& height, int& channels)
{
	fitsfile *fptr;
	int status = 0;
	fits_open_diskfile(&fptr, filename, READONLY, &status);
	if (status)
	{
		return status;
	}

	status = 0;
	int bitpix = 0, naxis = 0;
	long naxes[3] = { 0 };
	fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status);
	if (status)
	{
		status = 0;
		fits_close_file(fptr, &status);
		return status;
	}

	status = 0;
	double offset = 0;
	fits_read_key(fptr, TDOUBLE, "BZERO", &offset, NULL, &status);
	if (!status) {
		if (bitpix == SHORT_IMG && offset != 0.0) {
			bitpix = USHORT_IMG;
		}
		else if (bitpix == LONG_IMG && offset != 0.0)
			bitpix = ULONG_IMG;
	}
	else {
		/* but some software just put unsigned 16-bit data in the file
		* and don't set the BZERO keyword... */
		if (status == KEY_NO_EXIST && bitpix == SHORT_IMG)
			bitpix = USHORT_IMG;
	}

	int orig_bitpix = bitpix;

	int rx = naxes[0];
	int ry = naxes[1];
	int nbdata = rx * ry;

	

	if (naxis == 3 && naxes[2] != 3) {
		cout << "Unsupported FITS image with " << naxes[2] << " channels.\n";
		status = 0;
		fits_close_file(fptr, &status);
		return -1;
	}

	if (naxis == 2 && naxes[2] == 0) {
		naxes[2] = 1;
	}

	if (bitpix == LONGLONG_IMG) {
		cout<<"FITS images with 64 bits signed integer per pixel.channel are not supported.\n";
		status = 0;
		fits_close_file(fptr, &status);
		return -1;
	}

	width = rx;
	height = ry;
	channels = naxes[2];

	*img = new unsigned short[rx * ry * naxes[2]];
	if ((*img) <= 0)
	{
		cout << "new memery fail" << endl;
		fits_close_file(fptr, &status);
		return -1;
	}

	unsigned char *data8;
	double *pixels_double;
	uint32_t *pixels_long;
	long orig[3] = { 1L, 1L, 1L };
	// orig ^ gives the coordinate in each dimension of the first pixel to be read
	nbdata = naxes[0] * naxes[1] * naxes[2];

	fits_movabs_hdu(fptr, 1, 0, &status); // make sure reading primary HDU

	int zero = 0, datatype;

	switch (bitpix) {
	case BYTE_IMG:
		data8 = (unsigned char*)malloc(nbdata * sizeof(unsigned char));
		datatype = bitpix == BYTE_IMG ? TBYTE : TSBYTE;
		fits_read_pix(fptr, datatype, orig, nbdata, &zero,
			data8, &zero, &status);
		if (status) break;
		convert_data(bitpix, data8, *img, nbdata, FALSE);
		free(data8);
		break;
	case SHORT_IMG:
		fits_read_pix(fptr, TSHORT, orig, nbdata, &zero,
			*img, &zero, &status);
		if (status) break;
		//convert_data(bitpix, fit->data, fit->data, nbdata, FALSE);
		//fit->bitpix = USHORT_IMG;
		break;
	case USHORT_IMG:
		// siril 0.9 native, no conversion required
		fits_read_pix(fptr, TUSHORT, orig, nbdata, &zero,
			*img, &zero, &status);
		if (status == NUM_OVERFLOW) {
			// in case there are errors, we try short data
			status = 0;
			fits_read_pix(fptr, TSHORT, orig, nbdata, &zero, *img,
				&zero, &status);
			if (status)
				break;
			convert_data(SHORT_IMG, *img, *img, nbdata, FALSE);
			/*if (fit->lo)
				fit->lo += 32768;
			if (fit->hi)
				fit->hi += 32768;*/
			bitpix = USHORT_IMG;
		}
		break;
	case ULONG_IMG:		// 32-bit unsigned integer pixels
	case LONG_IMG:		// 32-bit signed integer pixels
		pixels_long = (uint32_t*)malloc(nbdata * sizeof(long));
		status = 0;
		datatype = bitpix == LONG_IMG ? TLONG : TULONG;
		fits_read_pix(fptr, datatype, orig, nbdata, &zero,
			pixels_long, &zero, &status);
		if (status) break;
		convert_data(bitpix, pixels_long, *img, nbdata, FALSE);
		free(pixels_long);
		bitpix = USHORT_IMG;
		break;
	case DOUBLE_IMG:	// 64-bit floating point pixels
	case FLOAT_IMG:		// 32-bit floating point pixels
		/*pixels_double = (double*)malloc(nbdata * sizeof(double));
		fits_read_pix(fptr, TDOUBLE, orig, nbdata, &zero,
			pixels_double, &zero, &status);
		if (status) break;
		convert_data(bitpix, pixels_double, *img, nbdata, fit->data_max > 1.0);
		free(pixels_double);
		fit->bitpix = USHORT_IMG;*/
		break;
	case LONGLONG_IMG:	// 64-bit integer pixels
	default:
		std::cout<<"FITS image format %d is not supported by Siril.\n";
		return -1;
	}
	status = 0;

	fits_close_file(fptr, &status);
	return status;
}

int create_catalog(std::vector<Edge*> edge, int width, int height, const char* filename)
{
	fitsfile *ofptr = NULL;
	int ncols = 4;
	int status = 0;
	char* ttype[] = { "X","Y","MAG_AUTO","FLUX", "LFLUX", "LBG" }; // { "X","Y","MAG_AUTO","FLUX", "LFLUX", "LBG" };
	char* tform[] = { "E","E","E","E","E","E" };
	char* tunit[] = { "pix","pix","unknown","unknown","unknown","unknown" };
	int nhdus, hdutype, nimgs;
	nimgs = 1;
	fits_create_file(&ofptr, filename, &status);
	fits_create_img(ofptr, 8, 0, NULL, &status);
	fits_write_key(ofptr, TSTRING, "SRCFN", (char*)filename, "Source image", &status);
	fits_create_tbl(ofptr, BINARY_TBL, edge.size(), ncols, ttype, tform,
		tunit, "SOURCES", &status);

	float *x = new float[edge.size()];
	int num = 0;
	for (auto val : edge)
	{
		x[num++] = val->x;
	}
	fits_write_col(ofptr, TFLOAT, 1, 1, 1, edge.size(), x, &status);
	num = 0;
	for (auto val : edge)
	{
		x[num++] = val->y;
	}
	fits_write_col(ofptr, TFLOAT, 2, 1, 1, edge.size(), x, &status);
	num = 0;
	for (auto val : edge)
	{
		x[num++] = val->val;
	}
	fits_write_col(ofptr, TFLOAT, 3, 1, 1, edge.size(), x, &status);
	num = 0;
	for (auto val : edge)
	{
		x[num++] = val->HFR;
	}
	fits_write_col(ofptr, TFLOAT, 4, 1, 1, edge.size(), x, &status);
	delete[] x;
	fits_modify_comment(ofptr, "TTYPE1", "X coordinate", &status);
	fits_modify_comment(ofptr, "TTYPE2", "Y coordinate", &status);
	fits_modify_comment(ofptr, "TTYPE3", "Sky background of source", &status);
	hdutype = IMAGE_HDU;
	fits_modify_comment(ofptr, "TTYPE4", "Flux of source", &status);

	fits_write_key(ofptr, TINT, "IMAGEW", &width, "Input image width", &status);

	fits_write_key(ofptr, TINT, "IMAGEH", &height, "Input image height", &status);

	fits_movabs_hdu(ofptr, 1, &hdutype, &status);
	assert(hdutype == IMAGE_HDU);
	fits_write_key(ofptr, TINT, "NEXTEND", &nimgs, "Number of extensions", &status);
	fits_close_file(ofptr, &status);
	return 0;
}


class NewStarDetection
{
private:
	double _hfd;
	std::vector<std::vector<double>> v_star_data;
public:
	cv::Point2f _BrightStarPosition;
	double _brightness;
	int _stars_num;
	std::vector<Edge*> deleteStars;
private:
	unsigned short* data = nullptr;
	unsigned short* dst = nullptr;
	int width;
	int height;
	unsigned short _light;
	std::vector<Edge*> v_edge;
private:
	inline unsigned short get_pixel(unsigned short* img, int y, int x)
	{
		return img[y* width + x];
	}
	inline void _dfs(int y, int x, std::vector<unsigned short>& v_pixels, std::vector<std::pair<int, int>>& pixel_coor)
	{
		std::list<std::pair<int, int>> list_pos;
		list_pos.push_back(std::pair<int, int>(y, x));
		while (!list_pos.empty())
		{
			int i = list_pos.back().first;
			int j = list_pos.back().second;
			list_pos.pop_back();
			/*unsigned short pixel = data[i * width + j];*/
			if (i >= 0 && j >= 0 && i < height && j < width && data[i * width + j] > _light)
			{
				v_pixels.push_back(data[i * width + j]);
				pixel_coor.push_back(std::pair<int, int>(i, j));
				data[i * width + j] = 0;

				list_pos.push_back(std::pair<int, int>(i - 1, j));
				list_pos.push_back(std::pair<int, int>(i, j - 1));
				list_pos.push_back(std::pair<int, int>(i, j + 1));
				list_pos.push_back(std::pair<int, int>(i + 1, j));
			}
		}
	}
	inline unsigned short _calc_light()
	{
		int len = width * height;
		double sum = std::accumulate(dst, dst + len, 0.0);
		double mean = sum / len; //均值

		double accum = 0.0;
		std::for_each(dst, dst + len, [mean, &accum](const double d) {
			accum += (d - mean)*(d - mean);
		});

		double stdev = sqrt(accum / (len - 1)); //方差

		_light = mean + 3 * stdev;
		return _light;
	}
	inline float halfFluxDiameter(std::vector<unsigned short>& v_pixel)
	{
		float totalFlux = 0;
		for (auto& val : v_pixel)
		{
			val = val - _light;
			totalFlux += val;
		}

		std::vector<unsigned short> v_stroe = v_pixel;

		float area = 0;
		float halfFlux = 0.5 * totalFlux;
		std::sort(v_stroe.begin(), v_stroe.end(), [](auto a, auto b) {
			return a > b;
		});

		//std::cout << "totalFlux : " << totalFlux << '\n';

		float nextAccum = 0., accum = 0;
		float last_val;
		for (auto val : v_stroe)
		{
			nextAccum = accum + val;
			last_val = val;
			if (nextAccum > halfFlux)
			{
				break;
			}
			accum += val;
			area += 1;
		}

		float deficit = halfFlux - accum;
		area += deficit / last_val;

		//std::cout << "area : " << area << '\n';

		float hfd = 2. *  sqrt(area / 3.141592653);

		return hfd;
	}
	inline std::pair<float, float> center_of_mass(std::vector<std::pair<int, int>>& pixel_coor, std::vector<unsigned short>& v_pixel)
	{
		float sum = 0, sumX = 0, sumY = 0;
		for (int i = 0, len = v_pixel.size(); i < len; i++)
		{
			sumX += pixel_coor[i].second * v_pixel[i];
			sumY += pixel_coor[i].first * v_pixel[i];
			sum += v_pixel[i];
		}

		float x = sumX / sum;
		float y = sumY / sum;
		return std::pair<float, float>(y, x);
	}
	inline int islandPerimeter(const std::vector<std::pair<int, int>>& pixel_coor)
	{
		int x0 = 100000, y0 = 100000;
		int x1 = 0, y1 = 0;
		for (auto it : pixel_coor)
		{
			int x = it.second;
			int y = it.first;

			x0 = x0 > x ? x : x0;
			x1 = x1 < x ? x : x1;
			y0 = y0 > y ? y : y0;
			y1 = y1 < y ? y : y1;
		}

		int count = 0;
		for (int i = y0; i <= y1; i++)
		{
			for (int j = x0; j <= x1; j++)
			{
				if (get_pixel(dst, i, j) > _light)
				{
					count += 4;
					if (i > y0 && get_pixel(dst, i - 1, j) > _light)
						count -= 2;
					if (j > x0 && get_pixel(dst, i, j - 1) > _light)
						count -= 2;
				}
			}
		}
		return count;
	}
	inline float compute_radius(const std::vector<std::pair<int, int>>& pixel_coor, int x, int y)
	{
		int x_num = 0, y_num = 0;
		for (auto it : pixel_coor)
		{
			if (it.first == y)
			{
				++y_num;
			}
			else if (it.second == x)
			{
				++x_num;
			}
		}

		if (x_num * 1. / y_num < 0.5 || x_num * 1. / y_num > 2.)
		{
			return 0;
		}
		else
		{
			return (x_num + y_num) * 1. / 4.;
		}
	}
	inline void getBoundingRectangle(const std::vector<std::pair<int, int>>& cloud, std::pair<int, int>& minXY, std::pair<int, int>& maxXY)
	{
		int minX = INT_MAX;
		int maxX = INT_MIN;
		int minY = INT_MAX;
		int maxY = INT_MIN;

		for (auto& pt : cloud)
		{
			int x = pt.second;
			int y = pt.first;

			// check X coordinate
			if (x < minX)
				minX = x;
			if (x > maxX)
				maxX = x;

			// check Y coordinate
			if (y < minY)
				minY = y;
			if (y > maxY)
				maxY = y;
		}
		if (minX > maxX) // if no point appeared to set either minX or maxX
			std::cout << "List of points can not be empty." << '\n';
		minXY = std::pair<int, int>(minY, minX);
		maxXY = std::pair<int, int>(maxY, maxX);
	}
	inline bool isCircle(const std::vector<std::pair<int, int>>& edgePoints)
	{
		float minAcceptableDistortion = 0.5f;
		float relativeDistortionLimit = 0.03f;

		float radius = 0.;

		if (edgePoints.size() < 8)
		{
			radius = 0;
			return false;
		}
		std::pair<int, int> minXY, maxXY;
		getBoundingRectangle(edgePoints, minXY, maxXY);

		std::pair<float, float> cloudSize(maxXY.first - minXY.first, maxXY.second - minXY.second);
		std::pair<float, float> t(cloudSize.first / 2, cloudSize.second / 2);

		std::pair<float, float> center(minXY.first + t.first, minXY.second + t.second);
		radius = ((float)cloudSize.second + cloudSize.first) / 4;

		float meanDistance = 0;
		for (int i = 0, n = edgePoints.size(); i < n; i++)
		{
			meanDistance += (float)abs(sqrt((center.first - edgePoints[i].first) * (center.first - edgePoints[i].first) +
				(center.second - edgePoints[i].second) * (center.second - edgePoints[i].second)) - radius);
		}

		meanDistance /= edgePoints.size();
		float maxDitance = std::max<float>(minAcceptableDistortion,
			((float)cloudSize.first + cloudSize.second) / 2 * relativeDistortionLimit);

		return (meanDistance <= maxDitance);
	}
	inline void detect_edge(const bool& bRun)
	{

		for (int i = 1; i < height - 1; i++)
		{
			for (int j = 1; j < width - 1; j++)
			{
				if (!bRun)
				{
					return;
				}
				if (get_pixel(data, i, j) > _light)
				{
					bool flag = true;
					//flag &= get_pixel(data, i - 1, j - 1) > _light;
					flag &= get_pixel(data, i - 1, j) > _light;
					//flag &= get_pixel(data, i - 1, j + 1) > _light;
					flag &= get_pixel(data, i, j - 1) > _light;
					flag &= get_pixel(data, i, j + 1) > _light;
					//flag &= get_pixel(data, i - 1, j + 1) > _light;
					flag &= get_pixel(data, i + 1, j) > _light;
					//flag &= get_pixel(data, i + 1, j + 1) > _light;

					if (flag)
					{
						std::vector<unsigned short> v_pixel;
						std::vector<std::pair<int, int>> pixel_coor;

						//检测星点区域
						_dfs(i, j, v_pixel, pixel_coor);
						//计算质心
						auto com = center_of_mass(pixel_coor, v_pixel);
						//计算星点半径
						float radius = compute_radius(pixel_coor, com.second, com.first);
						/*bool res = isCircle(pixel_coor);
						std::cout << res << '\n';*/

						float afa = 0; //afa计算
									   //					if (radius > 0)
						{
							float c = islandPerimeter(pixel_coor);
							float s = pixel_coor.size();
							afa = 4 * 3.141592653 * s / (c*c); //afa计算
						}

						if (radius > 0 && afa >= 0.1)
						{
							//计算星点 HFD
							float hfd = halfFluxDiameter(v_pixel);
							Edge * edge = new Edge();
							edge->x = com.second;
							edge->y = com.first;
							edge->HFR = hfd;
							edge->width = radius * 2;
							edge->val = dst[i * width + j];
							v_edge.push_back(edge);
						}
					}
				}
			}
		}
	}
	inline void setImage(const unsigned short * img, int width, int height, const bool& bRun)
	{
		v_edge.clear();
		if (dst != nullptr)
		{
			delete[] dst;
			dst = nullptr;
		}
		if (data != nullptr)
		{
			delete[] data;
			data = nullptr;
		}

		this->width = width;
		this->height = height;
		int len = width * height;
		data = new unsigned short[len];
		dst = new unsigned short[len];

		memcpy(data, img, sizeof(unsigned short) * len);
		cv::Mat src(height, width, CV_16UC1, (unsigned char*)data);
		cv::GaussianBlur(src, src, cv::Size(3, 3), 0.5, 0.5);
		memcpy(data, src.data, sizeof(unsigned short) * len);
		memcpy(dst, src.data, sizeof(unsigned short) * len);
		//计算检测星点阈值
		_calc_light();
		//星点检测
		detect_edge(bRun);
		std::sort(v_edge.begin(), v_edge.end(), [](auto a, auto b) {
			return a->HFR > b->HFR;
		});
	}
public:
	NewStarDetection();
	~NewStarDetection();
	void setmulsexImageData(unsigned short* image, int width, int height, float crop);
	double getAverageHFD();
	void clearPoint()
	{
		_BrightStarPosition.x = -1;
		_BrightStarPosition.y = -1;
	}
	std::vector<std::vector<double>> detectStarPoint(unsigned short* image, std::string bayer, int width, int height, bool &bRun);
	inline std::vector<Edge*> get_v_edge() { return v_edge; }
};


NewStarDetection::NewStarDetection()
{
	dst = nullptr;
	data = nullptr;
}

NewStarDetection::~NewStarDetection()
{
	if (data != nullptr)
	{
		delete[] data;
		data = nullptr;
	}
	if (dst != nullptr)
	{
		delete[] dst;
		dst = nullptr;
	}
}

void NewStarDetection::setmulsexImageData(unsigned short *image, int width, int height, float crop)
{
	int _width = static_cast<int>((1 - crop * 2) * width);
	int _height = static_cast<int>((1 - crop * 2) * height);
	int _x = crop * width;
	int _y = crop * height;
	_x = _x / 2 * 2;
	_y = _y / 2 * 2;
	bool brun = true;

	setImage(image, width, height, brun);
	std::cout << "v_edge num : " << v_edge.size() << std::endl;
	// std::vector<Edge*> p1;
	// std::vector<Edge*> p = v_edge;
	// for (unsigned int i = 0; i<p.size(); i++)
	// {
		// if (p[i]->HFR >= 1.5 && p[i]->HFR <= 40 && p[i]->x > _x&&p[i]->x < _x + _width&&p[i]->y > _y&&p[i]->y < _y + _height && p[i]->val < 60000)
		// {
			// p1.push_back(p[i]);
		// }
	// }
	// p.clear();
	// p = p1;
	// if (p.empty())
	// {
		// _hfd = 0;
		// _BrightStarPosition.x = -1;
		// _BrightStarPosition.y = -1;
		// _stars_num = 0;
	// }
	// else
	// {
		// //多星对焦，选择多颗星，排除前面10%的星，排除后面20%的星
		// int s_1 = 0.1 * p.size();
		// int s_2 = 0.8 * p.size();

		// if (s_2 - s_1 <= 0)
		// {
			// s_1 = 0;
			// s_2 = p.size();
		// }
		// std::cout << "---------------------------------" << '\n';
		// std::cout << "stars : " << p.size() << '\n';
		// std::cout << "mul stars : " << s_2 - s_1 << '\n';
		// std::cout << "---------------------------------" << '\n';
		// float sum_hfd = 0.;
		// for (int i = s_1; i < s_2; i++)
		// {
			// sum_hfd += p[i]->HFR;
		// }
		// _hfd = sum_hfd / (s_2 - s_1);
		// _stars_num = s_2 - s_1;
		// //取一个hfd最接近平均值的 星点用于显示

		// std::vector<std::pair<int, float> > v_p_dif;
		// deleteStars.clear();
		// for (int i = s_1; i < s_2; i++)
		// {
			// v_p_dif.push_back(std::pair<int, float>(i, fabs(p[i]->HFR - _hfd)));
			// deleteStars.push_back(p[i]);
		// }

		// std::sort(v_p_dif.begin(), v_p_dif.end(), [](std::pair<int, float> t1, std::pair<int, float> t2) {
			// return t1.second < t2.second;
		// });

		// int m_s = v_p_dif[0].first;
		// _BrightStarPosition.x = p[m_s]->x;
		// _BrightStarPosition.y = p[m_s]->y;
		// _brightness = p[m_s]->val;
		// std::cout << "display hfd : " << p[m_s]->HFR << '\n';
	// }
}

double NewStarDetection::getAverageHFD()
{
	return _hfd;
}

std::vector<std::vector<double> > NewStarDetection::detectStarPoint(unsigned short *image, std::string bayer, int width, int height, bool &bRun)
{
	std::vector<std::vector<double> > t;
	return t;
}

std::vector<Edge*> findStars(unsigned short * data, int width, int height)
{
	std::vector<Edge*> v_e;
	std::vector<std::pair<double, double>> v_p;
	cv::Mat img(height, width, CV_16UC1, data);
	double _mean = 0;
	double _std = 0;
	cv::Mat m1, s1;
	cv::meanStdDev(img, m1, s1);
	_mean = m1.at<double>(0, 0);
	_std = s1.at<double>(0, 0);

	double thresh = _mean + 3 * _std + 1;

	cv::Mat binarImg = img.clone();
	unsigned short* p = (unsigned short*)binarImg.data;

	int len = binarImg.cols * binarImg.rows;
	for (int i = 0; i < len; i++)
	{
		if (p[i] > thresh)
			p[i] = 65535;
		else
			p[i] = 0;
	}

	binarImg.convertTo(binarImg, 1. / 255);

	cv::Mat element1 = cv::getStructuringElement(cv::MorphShapes::MORPH_ELLIPSE, cv::Size(3, 3));
	cv::erode(binarImg, binarImg, element1);
	cv::Mat element2 = cv::getStructuringElement(cv::MorphShapes::MORPH_ELLIPSE, cv::Size(5, 5));
	cv::dilate(binarImg, binarImg, element2);

	unsigned char *p1 = binarImg.data;
	unsigned short *p2 = data;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if (p1[i * width + j] > 10)
			{
				std::list<std::vector<int>> lvi;
				p1[i * width + j] = 0;
				std::vector<int> v;
				v.push_back(i);
				v.push_back(j);
				v.push_back(p2[i * width + j]);
				lvi.push_back(v);

				double x = 0, y = 0, x_sum = 0, y_sum = 0;
				double sum = 0;

				while (!lvi.empty())
				{
					x_sum += lvi.front()[1] * lvi.front()[2];
					y_sum += lvi.front()[0] * lvi.front()[2];
					sum += lvi.front()[2];
					int kx = lvi.front()[1];
					int ky = lvi.front()[0];
					lvi.pop_front();
					//p1[ky * width + kx] = 0;

					if (kx - 1 > 0 && p1[(kx - 1) + ky * width] > 10)
					{
						v.clear();
						v.push_back(ky);
						v.push_back(kx - 1);
						v.push_back(p2[ky * width + kx - 1]);
						p1[ky * width + kx - 1] = 0;
						lvi.push_back(v);
					}

					if (kx + 1 < width && p1[(kx + 1) + ky * width] > 10)
					{
						v.clear();
						v.push_back(ky);
						v.push_back(kx + 1);
						v.push_back(p2[ky * width + kx + 1]);
						p1[ky * width + kx + 1] = 0;
						lvi.push_back(v);
					}

					if (ky - 1 > 0 && p1[kx + (ky - 1) * width] > 10)
					{
						v.clear();
						v.push_back(ky - 1);
						v.push_back(kx);
						v.push_back(p2[(ky - 1) * width + kx]);
						p1[(ky - 1) * width + kx] = 0;
						lvi.push_back(v);
					}

					if (ky + 1 < height && p1[kx + (ky + 1) * width] > 10)
					{
						v.clear();
						v.push_back(ky + 1);
						v.push_back(kx);
						v.push_back(p2[(ky + 1) * width + kx]);
						p1[(ky + 1) * width + kx] = 0;
						lvi.push_back(v);
					}
				}

				v_p.push_back(std::pair<double, double>(x_sum / sum, y_sum / sum));
				Edge e;
				e.x = x_sum / sum;
				e.y = y_sum / sum;
				e.val = img.at<unsigned short>(e.y, e.x);
				e.HFR = 3;
				v_e.push_back(&e);
			}
		}
	}

	//for (int i = 0; i < v_p.size(); i++)
	//{
	//	cv::circle(img, cv::Point(v_p[i].first, v_p[i].second), 10, cv::Scalar(65535), 1);
	//}

	//std::cout << v_p.size() << std::endl;

	//cv::imshow("b", img);
	//cv::imwrite("img.png", img);
	//cv::waitKey(0);
	std::cout << "findstars : " << v_e.size() << std::endl;
	return v_e;
}


int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		cout << "Incorrect input parameters" << '\n';
		return -1;
	}
	//char names[124] = "./pic/Stack_IC353_Light_10s_bin1_8.5C_gain100.fit";
	unsigned short* img = nullptr;
	int width = 0;
	int height = 0;
	int channels = 1;
	auto t1 = std::chrono::steady_clock::now();
	int status = readfits(argv[1], &img, width, height, channels);
	if (!status)
	{
		NewStarDetection star;
		//std::vector<Edge*> v_edge;
		if (channels == 3)
		{			
			//std::cout << channels << std::endl;
			star.setmulsexImageData(img + width * height, width, height, 0);
			//v_edge = findStars(img + width * height, width, height);
		}
		else
		{
			//std::cout << channels << std::endl;
			star.setmulsexImageData(img, width, height, 0);
			//v_edge = findStars(img, width, height);
		}
		
		//std::cout << channels << std::endl;
		std::vector<Edge*> v_edge = star.get_v_edge();
		
		create_catalog(v_edge, width, height, argv[2]);
	}
	
	if (img != nullptr)
	{
		delete[] img;
		img = nullptr;
	}
	auto t2 = std::chrono::steady_clock::now();

	std::chrono::duration<double, std::milli> mi = t2 - t1;
	std::cout << "detect stars sort time: " << mi.count() << "ms" << '\n';

	

	cout << "end" << endl;
	return 0;
}