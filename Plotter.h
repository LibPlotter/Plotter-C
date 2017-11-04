#pragma once
#include <stdlib.h>
#ifdef _WIN32
#include <glut.h>
#else 
#include <GL/glut.h>
#endif
#include <string>
#include <algorithm>

//Рисуем - поверхность?
bool is_surface = false;
//Параметры по умолчанию
//Текст
double def_text_size_big = 0.0075, def_text_size_normal = 0.005, def_text_size_small = 0.0025, def_text_size_very_small = 0.0015;
//Точки
double def_point_size_big = 5., def_point_size_normal = 2.5, def_point_size_small = 1.;
//Толщина линии
double def_line_thick = 2.0, def_line_normal = 1.0, def_line_thin = 0.5;
// Цвета по умолчанию
//Ось, текст
double color_axis[3], color_text[3];
//Динамический цвет
double color_array[3];
//Список цветов по умолчанию
int* color_table;
double *surface_colormap;
//Параметры изображения
//Отношение сторон
double aspect;
//Пределы изображения по осям
int add_x_canvas = 0, add_y_canvas = 0, add_z_canvas = 0;
//Ширина, высота, глубина
double diapason_width, diapason_height, diapason_depth;
//Начало и конец диапазонов отображения
double diapason_x_start, diapason_x_end, diapason_y_start, diapason_y_end, diapason_z_start, diapason_z_end;
//Здесь хранятся информация о осях для каждой поверхности - для корректной закраски
double **surface_start_end;
//Настройки вычисления - количество точек
int point_number_x, point_number_y;

//Экран
//Ширина и высота - передаются при инициализации
int window_width, window_height;
//Настройки вычисления
//Число наборов объектов, интервалов - передаются при инициализации в разных видах
int number_of_sets, number_of_intervals = 0;
//Шаг по осям
double step_x = 0.1, step_y = 0.1;
//Флаги освещения
bool lighting_activate = false, lighting_deactivate = false, lighting_status = false;
//Управление камерой
//Смещение по осям
int window_zoom = 0, x_shift = 0, y_shift = 0, z_shift = 0;
//Приближение/Удаление от точки
double factor = 1, new_factor = 1;
//Вращение
double x_degree = 0, y_degree = 0, z_degree = 0;
//Тип поверхности
//1 - точки
//2 - линии и точки
//3 - линии
//4 - поверхность
int surface_mode = 1;
//Тип графика линий
//1 - точечный
//2 - линии
//3 - линии с точками
int line_type = 2;
//Тип окрашивания поверхности
bool color_range = true;
//Параметры изображения
//Только оси
bool simplify3d = false;

//Строить засечки для оси Z ?
bool build_z_marks = false;
//Нужна ли легенда?
bool legend_need = false;
//Отображать сетку
bool net_need = true;

//Изображать все наборы
bool all_sets = false;
//Текущий отображаемый набор
int current_set = 0;
//Сглаживание
bool antialiasing = false;
//Показ подписей осей
bool axis_text_need = true;
//Миллиметровый режим
bool mm = false;

//Размер линий
//Ось
double axis_line_size = 2.0;
//Засечки
double axis_mark_line_size = 2.0;
//Сетка
double net_line_size = 0.5;
//График
double graphic_line_size = 3.0;
//Сеть поверхности
double surface_line_size = 5.0;
//Размер точек
//Осевые
double axis_point_size = 1.;
//График
double graphic_point_size = 5.;
//Поверхность
double surface_point_size = 5.;
//Толщина поверхности
double surface_triangle_volume = 0.5;
//Толщина линии текста
double text_size = 1.0;

//Надписи
std::string axis_ox_text = "", axis_oy_text = "", axis_oz_text = "";
//Размер надписей
double axis_text_size = 0.002, marks_text_size = 0.001;
//Относительное позиционирование
double x_x_position = 0, y_y_position = 0;
//Размер точек текста
double text_point_size = 0.1;
//Надписи легенды.
std::string* object_name;
//С какой позиции идёт легенда - от - 1 до 1
double zero_x = 0.66;

std::string color_name;
double * surface_tesseract;

//интерфейс функций
void draw_coded_character(int input, double coord_x, double coord_y, double size_coeff);
void background();
void draw_net();
void draw_axis();
void draw_axis_marks();
void draw();
void line_settings(double line_width, double color_r, double color_b, double  color_g);
void draw_line(double x_start, double y_start, double z_start, double x_end, double y_end, double z_end);
void draw_string(std::string input, double coord_x, double coord_y, double size_coeff);
int GetMagnitudeOfInteger(double x);
void draw_point(double x_pos, double y_pos, double z_pos);
void point_settings(double point_size, double color_r, double color_b, double color_g);
void keyboard(unsigned char key, GLint x, GLint y);
void initialization_bounds();
void set_bounds();
void scaling();
bool infinity_check(double evalue);
bool nan_check(double evalue);
void drawgraph_points_many_discontinuous(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number, double*** i_exclude_intervals, int i_exclude_pairs_number);
void drawgraph_points_many_perforated(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number, bool** i_perfocards);
void drawgraph_points_many_basic(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number);
void set_f_bounds();
void color_select(std::string color_name);
void factorization();
void draw_surface(int t);
void draw_graphic();
std::string color_map(int numer);
void drawgraph_surface_one_basic(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double i_diapason_y_start, double i_diapason_y_end, double(**i_surface_function)(double, double));
void drawgraph_surface_many_basic(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double i_diapason_y_start, double i_diapason_y_end, double(**i_surface_function)(double, double), int i_surface_number);
void surface_coloring(double axis_start, double axis_end, double current_axis);
void drawgraph_function_many_discontinuous(int i_window_width, int i_window_height, double  i_diapason_x_start, double i_diapason_x_end, double(**i_functions)(double), int i_functions_number, double ***i_exclude_intervals, int i_exclude_pairs_number);
void draw_string_3d(std::string input, double coord_x, double coord_y, double coord_z, double size_coeff);
void draw_string_stroke(std::string input);

class Point
{
public:
	Point()
	{
		valid = true;
		x = 0;
		y = 0;
		z = 0;
		for (int i = 0; i < 3; i++)
		{
			static_color[i] = 0;
			dynamic_color[i] = 0;
		}
	}
	bool valid;
	double x, y, z;
	double static_color[3];
	double dynamic_color[3];
};
void draw_triangle(Point first_point, Point second_point, Point third_point, double t_size);
class bound
{
public:
	Point left, right;
	bool eqv(int coord)
	{
		double temp;
		switch (coord)
		{
		case 0:
			temp = left.x - right.x;
			break;
		case 1:
			temp = left.y - right.y;
			break;
		case 2:
			temp = left.z - right.z;
			break;
		default:
			temp = 1;
		}
		temp = abs(temp);
		return temp < 0.00000001;
	}
};

bound** exclusion_intervals;
bool bounds_check(bound interval, double evalue, int coord);

class Graphic
{
public:
	Graphic(double(*new_fun)(double))
	{
		fun_y = new_fun;
		fun_z = new_fun;
		points = new Point[point_number_x + 1];
		Comp();
	}

	Graphic()
	{
		points = new Point[point_number_x + 1];
		fun_y = nullptr;
		fun_z = nullptr;

	}
	Point* points;
	double(*fun_y)(double);
	double(*fun_z)(double);

	//Расчёт точек поверхности
	void Comp()
	{
		for (int i = 0; i < point_number_x; i++)
		{
			bool skip = false;
			points[i].x = 0;
			points[i].y = 0;
			points[i].z = 0;
			double temp = diapason_x_start + step_x*i;
			for (int k = 0; k < number_of_intervals; k++)
				if (bounds_check(exclusion_intervals[current_set][k], temp, 0))
					skip = true;
			if (infinity_check(temp) || nan_check(temp))
				skip = true;

			if (skip == false)
			{
				points[i].x = temp;
				temp = fun_y(points[i].x);
				for (int k = 0; k < number_of_intervals; k++)
					if (bounds_check(exclusion_intervals[current_set][k], temp, 1))
						skip = true;
				if (infinity_check(temp) || nan_check(temp))
					skip = true;

			}
			if (skip == false)
			{
				points[i].y = temp;
				if (fun_z != nullptr)
				{
					temp = fun_z(points[i].x);
					for (int k = 0; k < number_of_intervals; k++)
						if (bounds_check(exclusion_intervals[current_set][k], temp, 2))
							skip = true;
					if (infinity_check(temp) || nan_check(temp))
						skip = true;
					points[i].z = temp;
				}
			}

			if (skip != false)	points[i].valid = false;
		}
	}
};

class Surface
{
public:
	//Поверхность
	double(*fun_surface)(double, double);
	Point **point;
	double surface_start, surface_end;
	Surface(double(*new_f)(double, double))
	{
		fun_surface = new_f;
		point = new Point*[point_number_x];
		for (int i = 0; i < point_number_x; i++)
			point[i] = new Point[point_number_y];

	}
	Surface()
	{
		point = new Point*[point_number_x];
		for (int i = 0; i < point_number_x; i++)
			point[i] = new Point[point_number_y]();
		fun_surface = nullptr;
	}
	void compute_surface()
	{
		surface_tesseract = new double[number_of_sets];
		double min_v = 1000000;
		double max_v = -1000000;
		for (int i = 0; i < point_number_x; i++)
			for (int j = 0; j < point_number_y; j++)
			{
				point[i][j].x = diapason_x_start + step_x*j;
				point[i][j].y = diapason_y_start + step_y*i;
				double temp = fun_surface(point[i][j].x, point[i][j].y);
				if (infinity_check(temp) || nan_check(temp))
				{
					point[i][j].z = 0;
					point[i][j].valid = true;
				}
				else
				{
					point[i][j].z = -1.*temp;
					if (temp > max_v)
						max_v = temp;
					if (temp < min_v)
						min_v = temp;
				}
			}
		surface_start = min_v;
		surface_end = max_v;
		diapason_z_start = std::min(min_v, diapason_z_start);
		diapason_z_start = floor(diapason_z_start);
		diapason_z_end = std::max(max_v, diapason_z_end);
		diapason_z_end = ceil(diapason_z_end);
	}
};

//Проверка: находится ли число между диапазонами.
//Возвращает TRUE если evalue лежит между left и right
bool bounds_check(bound interval, double evalue, int coord)
{
	if (interval.eqv(coord)) return false;
	switch (coord)
	{
	case 0:
		if (interval.left.x < evalue&&evalue < interval.right.x)
			return  true;
	case 1:
		if (interval.left.y < evalue&&evalue < interval.right.y)
			return  true;
	case 2:
		if (interval.left.z < evalue&&evalue < interval.right.z)
			return  true;
	default:
		return false;
	}
}
//Проверка на бесконечность
bool infinity_check(double evalue)
{
	if (std::isinf(evalue))
		return true;
	return false;
}
//Проверка nan
bool nan_check(double evalue)
{
	if (std::isnan(evalue))
		return true;
	return false;
}

double(**input_f)(double);
Point* points_array;
Graphic* function_array;
Surface* surf_array;

//Общая инициализация, передаются размеры окна
void initialization_general(int i_window_width, int i_window_height)
{
	window_width = i_window_width;
	window_height = i_window_height;

	color_select("Black");
	//Цвет текста
	color_text[0] = color_array[0];
	color_text[1] = color_array[1];
	color_text[2] = color_array[2];
	//Цвет осей
	color_select("Black");
	color_axis[0] = color_array[0];
	color_axis[1] = color_array[1];
	color_axis[2] = color_array[2];
	//Таблица цвета графиков по умолчанию
	color_table = new int[10];
	color_table[0] = 5;//Red
	color_table[1] = 6;//Lime
	color_table[2] = 15;//Navy
	color_table[3] = 16;//Magenta 
	color_table[4] = 11;//Cyan
	color_table[5] = 7;//Green
	color_table[6] = 0;//Gold
	color_table[7] = 1;//Yellow
	color_table[8] = 10;//Deep Sky Blue
	color_table[9] = 12;//Blue
	if (axis_ox_text == "")
		axis_ox_text = "X";
	if (axis_oy_text == "")
		axis_oy_text = "Y";
	if (axis_oz_text == "")
		axis_oz_text = "Z";
	if (object_name == nullptr)
	{
		object_name = new std::string[25];
		for (int i = 0; i < 25; i++)
		{
			object_name[i] = "Object N" + std::to_string(i + 1);

		}
	}
}

//Инициализация для функций
void initialization_function(double i_diapason_x_start, double i_diapason_x_end, double(**i_functions)(double), int i_functions_number)
{
	point_number_x = (i_diapason_x_end - i_diapason_x_start) / step_x + 1;
	number_of_sets = i_functions_number;
	function_array = new Graphic[number_of_sets];
	diapason_x_start = i_diapason_x_start;
	diapason_x_end = i_diapason_x_end;
	for (int i = 0; i < number_of_sets; i++)
	{
		current_set = i;
		function_array[i].fun_y = i_functions[i];
		function_array[i].Comp();
	}
	current_set = 0;
	//	initialization_bounds();
	set_f_bounds();
}

//Инициализация для поверхностей
void initialization_surface(double i_diapason_x_start, double i_diapason_x_end, double i_diapason_y_start, double i_diapason_y_end, double(**i_surface_functions)(double, double), int i_surface_number)
{
	point_number_x = (i_diapason_x_end - i_diapason_x_start) / step_x + 1;
	point_number_y = (i_diapason_y_end - i_diapason_y_start) / step_y + 1;
	number_of_sets = i_surface_number;
	diapason_x_start = i_diapason_x_start;
	diapason_x_end = i_diapason_x_end;
	x_x_position = i_diapason_x_end;
	diapason_y_start = i_diapason_y_start;
	diapason_y_end = i_diapason_y_end;
	y_y_position = i_diapason_y_end;
	surf_array = new Surface[number_of_sets];
	for (int i = 0; i < number_of_sets; i++)
	{
		current_set = i;
		surf_array[i].fun_surface = i_surface_functions[i];
		surf_array[i].compute_surface();
	}
}

//Инициализация интервала исключения
void initialization_exclusion(double*** i_exclude_intervals, int i_intervals_number, int i_function_number)
{
	number_of_intervals = i_intervals_number;
	number_of_sets = i_function_number;
	exclusion_intervals = new bound*[number_of_sets];
	for (int i = 0; i < number_of_sets; i++)
		exclusion_intervals[i] = new bound[number_of_intervals];
	for (int i = 0; i < number_of_intervals; i++)
	{
		for (int j = 0; j < number_of_sets; j++)
		{
			bound t;
			t.left.x = i_exclude_intervals[0][i][j];
			t.right.x = i_exclude_intervals[1][i][j];
			t.left.y = i_exclude_intervals[2][i][j];
			t.right.y = i_exclude_intervals[3][i][j];
			t.left.z = i_exclude_intervals[4][i][j];
			t.right.z = i_exclude_intervals[5][i][j];
			exclusion_intervals[j][i] = t;
		}
	}
}

//Инициализация наборов точек
void initialization_p_many(double ***i_points_xyz, int i_point_number_x, int i_axis_number, int i_sets_number)
{
	int arr_size = i_axis_number;
	point_number_x = i_point_number_x;
	number_of_sets = i_sets_number;
	function_array = new Graphic[number_of_sets];
	for (int i = 0; i < number_of_sets; i++)
		for (int j = 0; j < point_number_x; j++)
		{
			function_array[i].points[j].x = i_points_xyz[0][j][i];
			function_array[i].points[j].y = i_points_xyz[1][j][i];
			if (arr_size == 3)
				function_array[i].points[j].z = i_points_xyz[2][j][i];
			else
				function_array[i].points[j].z = 0;
		}
	initialization_bounds();
	set_bounds();
}

//Инициализация исключения точек
void initialization_p_exclude()
{
	for (int i = 0; i < number_of_sets; i++)
		for (int j = 0; j < point_number_x; j++)
			for (int k = 0; k < number_of_intervals; k++)
				//Лежит ли точка в одном интервале?
				if (bounds_check(exclusion_intervals[i][k], function_array[i].points[j].x, 0) || bounds_check(exclusion_intervals[i][k], function_array[i].points[j].y, 1) || bounds_check(exclusion_intervals[i][k], function_array[i].points[j].z, 2))
					function_array[i].points[j].valid = false;
}

//Задание размеров изображения графика для точеченого графика
void initialization_bounds()
{
	double temp = function_array[0].points[0].x;
	//min x
	for (int j = 0; j < number_of_sets; j++)
		for (int i = 0; i < point_number_x; i++)
		{
			if (function_array[j].points[i + 1].x < temp)
				temp = function_array[j].points[i + 1].x;
		}
	diapason_x_start = floor(temp);
	//max x
	for (int j = 0; j < number_of_sets; j++)
		for (int i = 0; i < point_number_x; i++)
		{
			if (function_array[j].points[i + 1].x > temp)
				temp = function_array[j].points[i + 1].x;
		}
	diapason_x_end = ceil(temp);
	x_x_position = diapason_x_end;
	temp = function_array[0].points[0].y;
	//min y
	for (int j = 0; j < number_of_sets; j++)
		for (int i = 0; i < point_number_x; i++)
		{
			if (function_array[j].points[i + 1].y < temp)
				temp = function_array[j].points[i + 1].y;
		}
	diapason_y_start = floor(temp);
	//max y
	for (int j = 0; j < number_of_sets; j++)
		for (int i = 0; i < point_number_x; i++)
		{
			if (function_array[j].points[i + 1].y > temp)
				temp = function_array[j].points[i + 1].y;
		}
	diapason_y_end = ceil(temp);
	y_y_position = diapason_y_end;
	//min z
	for (int j = 0; j < number_of_sets; j++)
		for (int i = 0; i < point_number_x; i++)
		{
			if (function_array[j].points[i + 1].z < temp)
				temp = function_array[j].points[i + 1].z;
		}
	diapason_y_start = floor(temp);
	set_bounds();
	//max z
	for (int j = 0; j < number_of_sets; j++)
		for (int i = 0; i < point_number_x; i++)
		{
			if (function_array[j].points[i + 1].z > temp)
				temp = function_array[j].points[i + 1].z;
		}
	diapason_y_end = ceil(temp);
}

//Задание размеров изображение после вычисления функции
void set_f_bounds()
{
	double min_q = std::min(diapason_x_start, diapason_x_end);
	diapason_x_start = min_q;
	diapason_y_start = min_q;
	double max_q = std::max(diapason_x_start, diapason_x_end);
	x_x_position = diapason_x_end;
	y_y_position = diapason_y_end;
	diapason_x_end = max_q;
	diapason_y_end = max_q;
	if (y_y_position > diapason_y_end) y_y_position = diapason_y_end;
	if (x_x_position > diapason_x_end) x_x_position = diapason_x_end;
}

//Задание размеров изображение после вычисления
void set_bounds()
{
	double min_q = std::min(std::min(diapason_x_start, diapason_x_end), std::min(diapason_y_start, diapason_y_end));
	diapason_x_start = min_q;
	diapason_y_start = min_q;
	double max_q = std::max(std::max(diapason_x_start, diapason_x_end), std::max(diapason_y_start, diapason_y_end));
	x_x_position = diapason_x_end;
	y_y_position = diapason_y_end;
	diapason_x_end = max_q;
	diapason_y_end = max_q;
	if (y_y_position > diapason_y_end) y_y_position = diapason_y_end;
	if (x_x_position > diapason_x_end) x_x_position = diapason_x_end;
}

//Передача информации о исключённых точках
void initialization_perfocard(int i_point_number_x, int i_set_number, bool** i_perfocards)
{
	for (int i = 0; i < i_set_number; i++)
		for (int j = 0; j < i_point_number_x; j++)
			function_array[i].points[j].valid = i_perfocards[i][j];
}

//Начальные настройки графики
void graphic_initialization()
{
	double w = window_width;
	double h = window_height;
	aspect = w / h;
	factorization();
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDepthMask(GL_TRUE);
	glClearColor(1.0, 1.0, 1.0, 0.0);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	if (is_surface == true)
	{
		glEnable(GL_POLYGON_SMOOTH);
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glShadeModel(GL_POLYGON_SMOOTH);
	}

	background();
}

//Легенда
void draw_legend()
{
	int from_set, to_set;
	if (all_sets == false)
	{
		from_set = current_set;
		to_set = current_set + 1;
	}
	else
	{
		from_set = 0;
		to_set = number_of_sets;
	}
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	color_select("Light Grey");
	line_settings(5., color_array[0], color_array[1], color_array[2]);
	double one = 1;
	double numbers = 0.06*(to_set - from_set + 1);
	double zero_y = one - numbers;
	glBegin(GL_QUADS);//from - 1, -1 to 1, 1
	glVertex3f(zero_x, zero_y, 0.);//0, 0	0, 0	->	1, 0
	glVertex3f(one, zero_y, 0.);//!1, 0 / \ |
	glVertex3f(one, one, 0.);// 1, 1 | \ /
	glVertex3f(zero_x, one, 0.);//0, 1	0, 1 < --	1, 1
	glEnd();
	for (int s = from_set; s < to_set; s++)
	{
		int temp = color_table[s];
		color_select(color_map(temp));
		line_settings(1., color_array[0], color_array[1], color_array[2]);
		glPushMatrix();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glBegin(GL_QUADS);//from - 1, -1 to 1, 1
		glVertex3f(zero_x + 0.01, one - 0.015 - 0.05*(s - from_set), 0.);//0, 0	0, 0	->	1, 0
		glVertex3f(zero_x + 0.04, one - 0.015 - 0.05*(s - from_set), 0.);//1, 0 / \ |
		glVertex3f(zero_x + 0.04, one - 0.055 - 0.05*(s - from_set), 0.);//1, 1 | \ /
		glVertex3f(zero_x + 0.01, one - 0.055 - 0.05*(s - from_set), 0.);//0, 1	0, 1 < --	1, 1
		glEnd();
		glTranslatef(zero_x + 0.05, one - 0.05*(s - from_set + 1), 0.);
		glScalef(0.0003, 0.0003, 0.0003);
		color_select("Black");
		line_settings(1., color_array[0], color_array[1], color_array[2]);
		draw_string_stroke(object_name[s]);
		glPopMatrix();
	}
	glPopMatrix();
}

//отрисовка кадра
void display()
{
	//Применение вращения, смешения и масштабирования
	factorization();
	//Фон
	background();
	//Сеть, оси и засечки. 
	//В упрощённом режиме только оси.
	if (simplify3d == false)
	{
		if (net_need == true)
			draw_net();
		draw_axis();
		draw_axis_marks();
	}
	else
		draw_axis();

	//отрисовка
	draw();
	if (legend_need == true)
		draw_legend();
	//вывод
	glutSwapBuffers();
}
//Вращение, смещение, факторизация
void factorization()
{
	double w = window_width;
	double h = window_height;
	aspect = w / h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, window_width, window_height);
	if (aspect > 1)
		glOrtho(double(diapason_x_start)*aspect, double(diapason_x_end)*aspect, double(diapason_y_start), double(diapason_y_end), -100., 100.);
	else
		glOrtho(double(diapason_x_start), double(diapason_x_end), double(diapason_y_start) / aspect, double(diapason_y_end) / aspect, -100., 100.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(1., 1., -1.);
	glRotatef(x_degree, 1., 0., 0.);
	glRotatef(y_degree, 0., 1., 0.);
	glRotatef(z_degree, 0., 0., 1.);
	glTranslatef(-x_shift / new_factor, -y_shift / new_factor, -z_shift / new_factor);
	glScalef(1 / new_factor, 1 / new_factor, 1 / new_factor);
}


//Фон
void background()
{
	glClear(GL_COLOR_BUFFER_BIT + GL_DEPTH_BUFFER_BIT);
}

//Сетка
void draw_net()
{
	//Полное смещение в плоскости
	double factor_x = 1;
	double factor_y = 1;
	if (aspect > 0)
		factor_x = aspect;
	else
		factor_y = aspect;
	int full_shift = abs(x_shift) + abs(y_shift);
	//Горизонтальная линия
	int end_point = (diapason_y_end - diapason_y_start + add_y_canvas + full_shift) * 2 / factor_y*factor_x;
	int start_point = -end_point;
	color_select("Dark Grey");
	line_settings(net_line_size, color_array[0], color_array[1], color_array[2]);
	for (int i = start_point; i < end_point; i++)
	{
		draw_line(double(start_point), double(i), 0., double(end_point), double(i), 0.);
		if (mm == true)
			for (int j = 0; j < 10; j++)
			{
				draw_line(double(start_point), double(i) + 0.1*double(j), 0., double(end_point), double(i) + 0.1*double(j), 0.);
			}
	}
	//Вертикальные линии
	end_point = (diapason_x_end - diapason_x_start + add_x_canvas + full_shift) * 2 / factor_y*factor_x;
	start_point = -end_point;

	for (int i = start_point; i < end_point; i++)
	{
		draw_line(double(i), double(start_point), 0., double(i), double(end_point), 0.);
		if (mm == true)
			for (int j = 1; j < 10; j++)
				draw_line(double(i) + 0.1*double(j), double(start_point), 0., double(i) + 0.1*double(j), double(end_point), 0.);
	}
}

//Отрисовка осей
void draw_axis()
{
	double full_shift = abs(x_shift) + abs(y_shift);
	double factor_x = 1;
	double factor_y = 1;
	double x_p = (abs(diapason_y_start) + abs(diapason_y_end) + (double)add_y_canvas + full_shift) / factor_y*factor_x;
	double y_p = (abs(diapason_x_start) + abs(diapason_x_end) + (double)add_x_canvas + full_shift) / factor_y*factor_x;
	double z_p = abs(diapason_z_start) + abs(diapason_z_end) + (double)add_z_canvas;
	line_settings(axis_line_size, color_axis[0], color_axis[1], color_axis[2]);
	draw_line(0., -x_p, 0., 0., x_p, 0.);
	draw_line(-y_p, 0., 0., y_p, 0., 0.);
	draw_line(0., 0., -5.*z_p, 0., 0., 5.*z_p);
	if (axis_text_need)
	{
		point_settings(axis_point_size / factor, color_axis[0], color_axis[1], color_axis[2]);
		line_settings(0.1, color_text[0], color_text[1], color_text[2]);
		draw_string(axis_ox_text, (diapason_x_end)*factor + 0.5, -0.5, axis_text_size);
		draw_string(axis_oy_text, -0.5, (diapason_y_end)*factor - 0.5, axis_text_size);
	}
}

//Отрисовка засечек
void draw_axis_marks()
{
	std::string minus;
	int net_points = 50;
	double factor_x = 1;
	double factor_y = 1;

	if (aspect > 0)
		factor_x = aspect;
	else
		factor_y = aspect;

	int full_shift = abs(x_shift) + abs(y_shift);
	//Порядок числа - сколько символов надо
	int string_size = GetMagnitudeOfInteger((double)net_points);


	int end_point = (diapason_y_end - diapason_y_start + add_y_canvas + full_shift) / factor_y*factor_x;
	int start_point = -end_point;
	for (int i = start_point; i < end_point; i++)
	{
		line_settings(axis_mark_line_size, color_axis[0], color_axis[1], color_axis[2]);
		draw_line(double(i), 0.2, 0., double(i), -0.2, 0.);


		point_settings(axis_point_size / factor, color_axis[0], color_axis[1], color_axis[2]);
		draw_point(double(i), 0., 0.);

		line_settings(text_size, color_text[0], color_text[1], color_text[2]);
		minus = std::to_string(i);
		draw_string(minus, double(i), 0.1, marks_text_size);
	}

	end_point = (diapason_x_end - diapason_x_start + add_x_canvas + full_shift) / factor_y*factor_x;
	start_point = -end_point;
	for (int i = start_point; i < end_point; i++)
	{
		line_settings(axis_mark_line_size, color_axis[0], color_axis[1], color_axis[2]);
		draw_line(0.2, double(i), 0., -0.2, double(i), 0.);


		point_settings(axis_point_size / factor, color_axis[0], color_axis[1], color_axis[2]);
		draw_point(0., double(i), 0.);

		line_settings(text_size, color_text[0], color_text[1], color_text[2]);
		minus = std::to_string(i);
		if (i != 0)
			draw_string(minus, 0.3, double(i), marks_text_size);
	}
	if (build_z_marks == true)
	{
		end_point = (diapason_z_end - diapason_z_start + abs(z_shift));
		start_point = -end_point;
		for (int i = start_point; i < end_point; i++)
		{
			line_settings(axis_mark_line_size, color_axis[0], color_axis[1], color_axis[2]);
			draw_line(-0.2, 0., double(i), 0.2, 0., double(i));

			color_select("Maroon");
			point_settings(axis_point_size / factor, color_axis[0], color_axis[1], color_axis[2]);
			draw_point(0., 0., double(i));

			line_settings(text_size, color_axis[0], color_axis[1], color_axis[2]);
			minus = std::to_string(i);
			if (i != 0)
			{
				draw_string_3d(minus, 0.2, 0.1, double(i), marks_text_size);
			}
		}
	}
}

//Отрисовка графика
void draw()
{
	if (is_surface == true)
	{
		//Точки
		if (surface_mode == 1 || surface_mode == 2)
			draw_surface(1);
		//Линии
		if (surface_mode == 2 || surface_mode == 3 || surface_mode == 4)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			draw_surface(2);
		}
		//Заливка
		if (surface_mode == 4)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			draw_surface(3);
		}
	}
	else
		draw_graphic();
}

void draw_surface(int t)
{
	int from_set, to_set;
	double t_size = surface_line_size;
	//Выбор наборов для отрисовки
	if (all_sets == false)
	{
		from_set = current_set;
		to_set = current_set + 1;
	}
	else
	{
		from_set = 0;
		to_set = number_of_sets;
	}

	if (t == 1)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		t_size = surface_point_size;
	}
	else if (t == 2)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		t_size = surface_line_size;
	}
	else if (t == 3)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		t_size = surface_triangle_volume;
	}
	t_size = t_size / factor;

	for (int s = from_set; s < to_set; s++)
	{
		int temp = color_table[s];
		color_select(color_map(temp));
		if (t == 2 || t == 3)
			line_settings(t_size, color_array[0], color_array[1], color_array[2]);
		else
			point_settings(t_size, color_array[0], color_array[1], color_array[2]);
		for (int i = 1; i < point_number_x; i++)
			for (int j = 1; j < point_number_y; j++)
			{
				//Проверка на скрытость точек
				if (surf_array[s].point[i][j].valid == true && surf_array[s].point[i - 1][j].valid == true && surf_array[s].point[i - 1][j - 1].valid == true)
					draw_triangle(surf_array[s].point[i][j], surf_array[s].point[i - 1][j], surf_array[s].point[i - 1][j - 1], t_size);
				if (surf_array[s].point[i][j].valid == true && surf_array[s].point[i][j - 1].valid == true && surf_array[s].point[i - 1][j - 1].valid == true)
					draw_triangle(surf_array[s].point[i][j], surf_array[s].point[i][j - 1], surf_array[s].point[i - 1][j - 1], t_size);
			}
	}

}

void draw_graphic()
{
	int from_set, to_set;
	if (all_sets == false)
	{
		from_set = current_set;
		to_set = current_set + 1;
	}
	else
	{
		from_set = 0;
		to_set = number_of_sets;
	}
	for (int s = from_set; s < to_set; s++)
	{
		int temp = color_table[s];
		color_select(color_map(temp));
		for (int i = 1; i < point_number_x; i++)
		{
			if (function_array[s].points[i - 1].valid == true && function_array[s].points[i].valid == true)
			{
				if (line_type != 2)
				{
					point_settings(graphic_point_size / factor, color_array[0], color_array[1], color_array[2]);
					draw_point(function_array[s].points[i - 1].x, function_array[s].points[i - 1].y, -function_array[s].points[i - 1].z);
				}
				if (line_type != 1)
				{
					line_settings(graphic_line_size, color_array[0], color_array[1], color_array[2]);
					draw_line(function_array[s].points[i - 1].x, function_array[s].points[i - 1].y, -function_array[s].points[i - 1].z, function_array[s].points[i].x, function_array[s].points[i].y, -function_array[s].points[i].z);
				}
			}
		}
		if (line_type != 2)
		{
			point_settings(graphic_point_size / factor, color_array[0], color_array[1], color_array[2]);
			draw_point(function_array[s].points[point_number_x - 1].x, function_array[s].points[point_number_x - 1].y, -function_array[s].points[point_number_x - 1].z);
		}
	}
}

//Параметры линии
void line_settings(double line_width, double color_r, double color_b, double  color_g)
{
	glLineWidth(line_width);
	glColor4f(color_r, color_b, color_g, 1.);
}

//Параметры точек
void point_settings(double point_size, double color_r, double color_b, double color_g)
{
	glPointSize(point_size);
	glColor4f(color_r, color_b, color_g, 1.);
}

//Отрисовка точек
void draw_point(double x_pos, double y_pos, double z_pos)
{
	glBegin(GL_POINTS);
	glVertex3f(x_pos, y_pos, z_pos);
	glEnd();
}

//Отрисовка линии
void draw_line(double x_start, double y_start, double z_start, double x_end, double y_end, double z_end)
{
	glBegin(GL_LINES);
	glVertex3f(x_start, y_start, z_start);
	glVertex3f(x_end, y_end, z_end);
	glEnd();
}

//Отрисовка треугольника
void draw_triangle(Point first_point, Point second_point, Point third_point, double t_size)
{
	glBegin(GL_TRIANGLES);
	glBegin(GL_TRIANGLES);
	if (all_sets == false)
	{
		surface_coloring(surf_array[current_set].surface_start, surf_array[current_set].surface_end, first_point.z);
		line_settings(t_size, color_array[0], color_array[1], color_array[2]);
	}
	glVertex3f(first_point.x, first_point.y, first_point.z);

	if (all_sets == false)
	{
		surface_coloring(surf_array[current_set].surface_start, surf_array[current_set].surface_end, second_point.z);
		line_settings(t_size, color_array[0], color_array[1], color_array[2]);
	}
	glVertex3f(second_point.x, second_point.y, second_point.z);

	if (all_sets == false)
	{
		surface_coloring(surf_array[current_set].surface_start, surf_array[current_set].surface_end, third_point.z);
		line_settings(t_size, color_array[0], color_array[1], color_array[2]);
	}
	glVertex3f(third_point.x, third_point.y, third_point.z);
	glEnd();
}

//Получение порядка числа
int GetMagnitudeOfInteger(double x)
{
	double temp = x;
	int y = 0;
	do
	{
		temp = temp / 10;
		y = y + 1;
	} while (temp > 0);
	return y;
}

//Отрисовка текста - подготовка
void draw_text_preparation(double coord_x, double coord_y, double coord_z, double size_coeff)
{
	double w = window_width;
	double h = window_height;
	aspect = w / h;
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (aspect > 1)
		glOrtho(double(diapason_x_start)*aspect, double(diapason_x_end)*aspect, double(diapason_y_start), double(diapason_y_end), double(-100), double(100));
	else
		glOrtho(double(diapason_x_start), double(diapason_x_end), double(diapason_y_start) / aspect, double(diapason_y_end) / aspect, double(-100), double(100));
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(x_degree, 1., 0., 0.);
	glRotatef(y_degree, 0., 1., 0.);
	glRotatef(z_degree, 0., 0., 1.);
	glTranslatef((coord_x - x_shift) / new_factor, (coord_y - y_shift) / new_factor, (coord_z - z_shift) / new_factor);
	glScalef(size_coeff / new_factor, size_coeff / new_factor, size_coeff / new_factor);
}

double color_computing(double x)
{
	double temp = abs(x);
	do
	{
		temp = temp / 10;
	} while (temp > 1);
	return temp;
}

//Пишем текст
void draw_string_stroke(std::string input)
{
	glPointSize(text_point_size);
	for (int i = 0; i < input.size(); i++) glutStrokeCharacter(GLUT_STROKE_ROMAN, (int)input[i]);
}

void draw_string_3d(std::string input, double coord_x, double coord_y, double coord_z, double size_coeff)
{

	draw_text_preparation(coord_x, coord_y, coord_z, size_coeff);
	draw_string_stroke(input);
	glPopMatrix();
}
//Полный текст
void draw_string(std::string input, double coord_x, double coord_y, double size_coeff)
{
	draw_string_3d(input, coord_x, coord_y, 0.0, size_coeff);
}

//Изменение размера
void OnRedraw(GLint w, GLint h)
{
	window_width = w;
	window_height = h;

}

//Вызов отрисовки
void OpenGL()
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(window_width, window_height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("C++ Plotter");
	graphic_initialization();
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(OnRedraw);
	glutMainLoop();
}

//Получить цвет по номеру
std::string color_map(int numer)
{
	std::string color;
	switch (numer)
	{
	case 0:
	{
		color = "Gold";
		break;
	}
	case 1:
	{
		color = "Yellow";
		break;
	}
	case 2:
	{
		color = "Maroon";
		break;
	}
	case 3:
	{
		color = "Crimson";
		break;
	}
	case 4:
	{
		color = "Dark Red";
		break;
	}
	case 5:
	{
		color = "Red";
		break;
	}
	case 6:
	{
		color = "Lime";
		break;
	}
	case 7:
	{
		color = "Green";
		break;
	}
	case 8:
	{
		color = "Spring Green";
		break;
	}
	case 9:
	{
		color = "Dark Green";
		break;
	}
	case 10:
	{
		color = "Deep Sky Blue";
		break;
	}
	case 11:
	{
		color = "Cyan";
		break;
	}
	case 12:
	{
		color = "Blue";
		break;
	}
	case 13:
	{
		color = "Medium Blue";
		break;
	}
	case 14:
	{
		color = "Dark Blue";
		break;
	}
	case 15:
	{
		color = "Navy";
		break;
	}
	case 16:
	{
		color = "Magenta";
		break;
	}
	case 17:
	{
		color = "Black";
		break;
	}
	case 18:
	{
		color = "Silver";
		break;
	}
	case 19:
	{
		color = "Light Grey";
		break;
	}
	case 20:
	{
		color = "Dark Grey";
		break;
	}
	case 21:
	{
		color = "Grey";
		break;
	}
	default:color = "Red";
	}
	return color;
}

//Получение цвета по названию
void color_select(std::string color_name)
{
	color_array[0] = 0.;
	color_array[1] = 0.;
	color_array[2] = 0.;
	if (color_name == "Gold")
	{
		color_array[0] = 1;
		color_array[1] = 0.843;
	}
	if (color_name == "Yellow")
	{
		color_array[0] = 1;
		color_array[1] = 1;
	}
	if (color_name == "Maroon")
		color_array[0] = 0.502;
	if (color_name == "Crimson")
	{
		color_array[0] = 0.863;
		color_array[1] = 0.078;
		color_array[2] = 0.235;
	}
	if (color_name == "Dark Red")
		color_array[0] = 0.545;
	if (color_name == "Red")
		color_array[0] = 1;
	if (color_name == "Lime")
		color_array[1] = 1;
	if (color_name == "Green")
		color_array[1] = 0.502;
	if (color_name == "Spring Green")
	{
		color_array[1] = 1;
		color_array[2] = 0.498;
	}
	if (color_name == "Dark Green")
		color_array[1] = 0.392;
	if (color_name == "Deep Sky Blue")
	{
		color_array[1] = 0.749;
		color_array[2] = 1;
	}
	if (color_name == "Cyan")
	{
		color_array[1] = 1;
		color_array[2] = 1;
	}
	if (color_name == "Blue")
		color_array[2] = 1;
	if (color_name == "Medium Blue")
		color_array[2] = 0.804;
	if (color_name == "Dark Blue")
		color_array[2] = 0.545;
	if (color_name == "Navy")
		color_array[2] = 0.502;
	if (color_name == "Magenta")
	{
		color_array[0] = 1;
		color_array[2] = 1;
	}
	if (color_name == "Black")
	{
		color_array[0] = 0;
		color_array[1] = 0;
		color_array[2] = 0;
	}
	if (color_name == "Silver")
	{
		color_array[0] = 0.753;
		color_array[1] = 0.753;
		color_array[2] = 0.753;
	}
	if (color_name == "Light Grey")
	{
		color_array[0] = 0.827;
		color_array[1] = 0.827;
		color_array[2] = 0.827;
	}
	if (color_name == "Dark Grey")
	{
		color_array[0] = 0.663;
		color_array[1] = 0.663;
		color_array[2] = 0.663;
	}
	if (color_name == "Grey")
	{
		color_array[0] = 0.502;
		color_array[1] = 0.502;
		color_array[2] = 0.502;
	}
}

//Цвет поверхности
//Преобразуем точку y в точку между 0 и 1
//На основе этого вычисляем цвет
void surface_coloring(double axis_start, double axis_end, double current_axis)
{
	double length = axis_end - axis_start;
	double c_axis = abs(current_axis - axis_start) / length;
	if (c_axis <= current_axis && c_axis < 0.2)
	{
		color_select("Magenta");
		color_array[0] = 1 - (c_axis * 5);
	}
	else if (c_axis <= 0.2 && c_axis < 0.4)
	{
		color_select("Blue");
		color_array[1] = (c_axis - 0.2) * 5;
	}
	else if (c_axis <= 0.4 && c_axis < 0.6)
	{
		color_select("Cyan");
		color_array[2] = 1 - (c_axis - 0.4) * 5;
	}
	else if (c_axis <= 0.6 && c_axis < 0.8)
	{
		color_select("Green");
		color_array[0] = (c_axis - 0.6) * 5;
	}
	else if (c_axis <= 0.8 && c_axis < 1)
	{
		color_select("Yellow");
		color_array[1] = 1 - (c_axis - 0.8) * 5;
	}
	else
		color_select("Red");
}

//Конвертер отрисовки одной поверхности
void drawgraph_surface_one_basic(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double i_diapason_y_start, double i_diapason_y_end, double(**i_surface_function)(double, double))
{
	double(** c_surface_function)(double, double) = new(double(*[1])(double, double));
	//Конвертируем в массив функций
	c_surface_function = i_surface_function;
	drawgraph_surface_many_basic(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, i_diapason_y_start, i_diapason_y_end, c_surface_function, 1);
}


//Инициализация отрисовки многих поверхностей
void drawgraph_surface_many_basic(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double i_diapason_y_start, double i_diapason_y_end, double(**i_surface_function)(double, double), int i_surface_number)
{
	initialization_general(i_window_width, i_window_height);
	is_surface = true;
	initialization_surface(i_diapason_x_start, i_diapason_x_end, i_diapason_y_start, i_diapason_y_end, i_surface_function, i_surface_number);
	OpenGL();
}

//Инициализация отрисовки многих функций
void drawgraph_function_many_basic(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double(**i_functions)(double), int i_functions_number)
{
	initialization_general(i_window_width, i_window_height);
	initialization_function(i_diapason_x_start, i_diapason_x_end, i_functions, i_functions_number);
	OpenGL();
}
//Конвертер отрисовки одной функции
void drawgraph_function_one_basic(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double(**i_function)(double))
{
	double(**c_function)(double) = new(double(*[1])(double));
	c_function = i_function;
	drawgraph_function_many_basic(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, c_function, 1);
}

//Конвертер отрисовки одной функии с ограничениями
void drawgraph_function_one_discontinuous(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double(**i_function)(double), double** i_exclude_intervals, int i_exclude_pairs_number)
{
	double***c_exclude_intervals = new double**[6];
	for (int i = 0; i < 6; i++)
	{
		c_exclude_intervals[i] = new double*[i_exclude_pairs_number];
		for (int j = 0; j < i_exclude_pairs_number; j++)
		{
			c_exclude_intervals[i][j] = new double[i_exclude_pairs_number];
		}
	}

	double(**c_function)(double) = new(double(*[1])(double));
	c_function = i_function;
	for (int i = 0; i < i_exclude_pairs_number; i++)
		for (int j = 0; j < 6; j++)
			c_exclude_intervals[j][i][0] = i_exclude_intervals[j][i];
	drawgraph_function_many_discontinuous(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, c_function, 1, c_exclude_intervals, i_exclude_pairs_number);
}

//Инициализация отрисовки многих функций с ограничениями	
void drawgraph_function_many_discontinuous(int i_window_width, int i_window_height, double  i_diapason_x_start, double i_diapason_x_end, double(**i_functions)(double), int i_functions_number, double ***i_exclude_intervals, int i_exclude_pairs_number)
{
	initialization_general(i_window_width, i_window_height);
	initialization_exclusion(i_exclude_intervals, i_exclude_pairs_number, i_functions_number);
	initialization_function(i_diapason_x_start, i_diapason_x_end, i_functions, i_functions_number);
	OpenGL();
}
//Конвертер отрисовки одного точечного графика
void drawgraph_points_one_basic(int i_window_width, int i_window_height, double **i_points_xyz, int i_point_number, int i_axis_number)
{

	int arr_size = i_axis_number;
	double***c_points_xyz = new double**[arr_size];
	for (int i = 0; i < arr_size; i++)
	{
		c_points_xyz[i] = new double*[i_point_number];
		for (int j = 0; j < i_point_number; j++)
		{
			c_points_xyz[i][j] = new double[i_point_number];
		}
	}
	for (int i = 0; i < i_point_number; i++)
	{
		c_points_xyz[0][i][0] = i_points_xyz[0][i];
		c_points_xyz[1][i][0] = i_points_xyz[1][i];
		if (arr_size == 3)
			c_points_xyz[2][i][0] = i_points_xyz[2][i];

	}
	drawgraph_points_many_basic(i_window_width, i_window_height, c_points_xyz, i_point_number, i_axis_number, 1);
}


//Конвертер отрисовки одного точечного графика с точками исключения
void drawgraph_points_one_perforated(int i_window_width, int i_window_height, double ** i_points_xyz, int i_point_number, int i_axis_number, bool* i_perfocard)
{
	int arr_size = i_axis_number;
	double***c_points_xyz = new double**[arr_size];
	for (int i = 0; i < arr_size; i++)
	{
		c_points_xyz[i] = new double*[i_point_number];
		for (int j = 0; j < i_point_number; j++)
		{
			c_points_xyz[i][j] = new double[i_point_number];
		}
	}
	bool**c_perfocard = new bool*[1];
	for (int i = 0; i < 1; i++)
		c_perfocard[i] = new bool[i_point_number];
	for (int i = 0; i < i_point_number; i++)
	{
		c_points_xyz[0][i][0] = i_points_xyz[0][i];
		c_points_xyz[1][i][0] = i_points_xyz[1][i];
		if (arr_size == 3)
			c_points_xyz[2][i][0] = i_points_xyz[2][i];
		c_perfocard[0][i] = i_perfocard[i];
	}
	drawgraph_points_many_perforated(i_window_width, i_window_height, c_points_xyz, i_point_number, i_axis_number, 1, c_perfocard);
}

//Конвертер отрисовки одного точечного графика с интервалами исключения
void drawgraph_points_one_discontinuous(int i_window_width, int i_window_height, double** i_points_xyz, int i_point_number, int i_axis_number, double** i_exclude_intervals, int i_exclude_pairs_number)
{
	int arr_size = i_axis_number;
	double***c_points_xyz = new double**[arr_size];
	for (int i = 0; i < arr_size; i++)
	{
		c_points_xyz[i] = new double*[i_point_number];
		for (int j = 0; j < i_point_number; j++)
		{
			c_points_xyz[i][j] = new double[i_point_number];
		}
	}
	double***c_exclude_intervals = new double**[6];
	for (int i = 0; i < 6; i++)
	{
		c_exclude_intervals[i] = new double*[i_exclude_pairs_number];
		for (int j = 0; j < i_exclude_pairs_number; j++)
		{
			c_exclude_intervals[i][j] = new double[i_exclude_pairs_number];
		}
	}
	for (int i = 0; i < i_point_number; i++)
	{
		c_points_xyz[0][i][0] = i_points_xyz[0][i];
		c_points_xyz[1][i][0] = i_points_xyz[1][i];
		if (arr_size == 3)
			c_points_xyz[2][i][0] = i_points_xyz[2][i];
	}
	for (int i = 0; i < i_exclude_pairs_number; i++)
		for (int j = 0; j < 6; j++)
			c_exclude_intervals[j][i][0] = i_exclude_intervals[j][i];
	drawgraph_points_many_discontinuous(i_window_width, i_window_height, c_points_xyz, i_point_number, i_axis_number, 1, c_exclude_intervals, i_exclude_pairs_number);
}

//Инициализация отрисовки нескольких точечных графиков
void drawgraph_points_many_basic(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number)
{
	initialization_general(i_window_width, i_window_height);
	initialization_p_many(i_points_xyz, i_point_number, i_axis_number, i_set_number);
	OpenGL();
}

//Инициализация отрисовки нескольких точечных графиков
void drawgraph_points_many_perforated(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number, bool** i_perfocards)
{
	initialization_general(i_window_width, i_window_height);
	initialization_p_many(i_points_xyz, i_point_number, i_axis_number, i_set_number);
	initialization_perfocard(i_point_number, i_set_number, i_perfocards);
	OpenGL();
}

//Инициализация отрисовки нескольких точечных графиков
void drawgraph_points_many_discontinuous(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number, double*** i_exclude_intervals, int i_exclude_pairs_number)
{
	initialization_general(i_window_width, i_window_height);
	initialization_p_many(i_points_xyz, i_point_number, i_axis_number, i_set_number);
	initialization_exclusion(i_exclude_intervals, i_exclude_pairs_number, i_set_number);
	initialization_p_exclude();
	OpenGL();
}

void Draw(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double(**i_function)(double))
{
	drawgraph_function_one_basic(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, i_function);
}

void Draw(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double(**i_function)(double), double** i_exclude_intervals, int i_exclude_pairs_number)
{
	drawgraph_function_one_discontinuous(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, i_function, i_exclude_intervals, i_exclude_pairs_number);
}

void Draw(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double(**i_functions)(double), int i_functions_number)
{
	drawgraph_function_many_basic(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, i_functions, i_functions_number);
}

void Draw(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double(**i_functions)(double), int i_functions_number, double ***i_exclude_intervals, int i_exclude_pairs_number)
{
	drawgraph_function_many_discontinuous(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, i_functions, i_functions_number, i_exclude_intervals, i_exclude_pairs_number);
}

void Draw(int i_window_width, int i_window_height, double **i_points_xyz, int i_point_number, int i_axis_number)
{
	drawgraph_points_one_basic(i_window_width, i_window_height, i_points_xyz, i_point_number, i_axis_number);
}

void Draw(int i_window_width, int i_window_height, double ** i_points_xyz, int i_point_number, int i_axis_number, bool* i_perfocard)
{
	drawgraph_points_one_perforated(i_window_width, i_window_height, i_points_xyz, i_point_number, i_axis_number, i_perfocard);
}

void Draw(int i_window_width, int i_window_height, double** i_points_xyz, int i_point_number, int i_axis_number, double** i_exclude_intervals, int i_exclude_pairs_number)
{
	drawgraph_points_one_discontinuous(i_window_width, i_window_height, i_points_xyz, i_point_number, i_axis_number, i_exclude_intervals, i_exclude_pairs_number);
}

void Draw(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number)
{
	drawgraph_points_many_basic(i_window_width, i_window_height, i_points_xyz, i_point_number, i_axis_number, i_set_number);
}

void Draw(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number, bool** i_perfocards)
{
	drawgraph_points_many_perforated(i_window_width, i_window_height, i_points_xyz, i_point_number, i_axis_number, i_set_number, i_perfocards);
}

void Draw(int i_window_width, int i_window_height, double*** i_points_xyz, int i_point_number, int i_axis_number, int i_set_number, double*** i_exclude_intervals, int i_exclude_pairs_number)
{
	drawgraph_points_many_discontinuous(i_window_width, i_window_height, i_points_xyz, i_point_number, i_axis_number, i_set_number, i_exclude_intervals, i_exclude_pairs_number);
}

void Draw(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double i_diapason_y_start, double i_diapason_y_end, double(**i_surface_function)(double, double))
{
	drawgraph_surface_one_basic(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, i_diapason_y_start, i_diapason_y_end, i_surface_function);
}
void Draw(int i_window_width, int i_window_height, double i_diapason_x_start, double i_diapason_x_end, double i_diapason_y_start, double i_diapason_y_end, double(**i_surface_function)(double, double), int i_surface_number)
{
	drawgraph_surface_many_basic(i_window_width, i_window_height, i_diapason_x_start, i_diapason_x_end, i_diapason_y_start, i_diapason_y_end, i_surface_function, i_surface_number);
}


double compute_factor_coeff(double current_factor, double wanted_factor)
{
	return current_factor / wanted_factor;
}

//Замыкание круга
double degree_circling(double x)
{
	double y;
	if (x >= 360)
		y = x - 360;
	else
		y = x;
	return y;
}

//Обработчик клавиатуры. 
void keyboard(unsigned char key, GLint x, GLint y)
{
	switch (key)
	{
	case 'u':
	{
		if (antialiasing == false)
		{
			glEnable(GL_LINE_SMOOTH);
			glHint(GL_LINE_SMOOTH, GL_NICEST);
		}
		else

			glDisable(GL_LINE_SMOOTH);
		antialiasing = !antialiasing;

		break;
	}
	case 'o':
	{
		if (current_set < number_of_sets - 1)
			current_set = current_set + 1;
		else
			current_set = 0;
		break;
	}
	case 'i':
	{
		all_sets = !all_sets;
		break;
	}
	case't'://Переключение режима изображения
	{
		if (is_surface == false)
			if (line_type < 3)
				line_type = line_type + 1;
			else
				line_type = 1;
		else
			if (surface_mode < 4)
				surface_mode = surface_mode + 1;
			else
				surface_mode = 1;
		break;
	}
	case 'x':
	{
		x_degree = x_degree + 1;
		break;
	}
	case 'X':
	{
		x_degree = x_degree + 359;
		break;
	}
	case 'y':
	{
		y_degree = y_degree + 1;
		break;
	}
	case 'Y':
	{
		y_degree = y_degree + 359;
		break;
	}
	case 'z':
	{
		z_degree = z_degree + 1;
		break;
	}
	case 'Z':
	{
		z_degree = z_degree + 359;
		break;
	}
	case '-':
	{
		factor = new_factor;
		new_factor = factor + 0.01;
		add_x_canvas = add_x_canvas + 1;
		add_y_canvas = add_y_canvas + 1;
		break;
	}
	case '+':
	{
		if (factor > 0.05)
		{
			factor = new_factor;
			new_factor = factor - 0.01;
		}
		if (add_x_canvas > 0)
			add_x_canvas = add_x_canvas - 1;

		if (add_y_canvas > 0)
			add_y_canvas = add_y_canvas - 1;
		break;

	}
	case 'w':
	{
		y_shift = y_shift + 1;
		break;
	}
	case 's':
	{
		y_shift = y_shift - 1;
		break;
	}
	case 'a':
	{
		x_shift = x_shift - 1;
		break;
	}
	case 'd':
	{
		x_shift = x_shift + 1;
		break;
	}
	case 'q':
	{
		z_shift = z_shift + 1;
		break;
	}
	case 'e':
	{
		z_shift = z_shift - 1;
		break;
	}
	case 'l':
	{
		legend_need = !legend_need;
		break;
	}
	case 'n':
	{
		net_need = !net_need;
		break;
	}
	case 'm':
	{
		mm = !mm;
		break;
	}
	case 'A':
	{
		axis_text_need = !axis_text_need;
		break;
	}
	default:break;
	}
	//Обработка замыкания круга
	x_degree = degree_circling(x_degree);
	y_degree = degree_circling(y_degree);
	z_degree = degree_circling(z_degree);
	display();
}

//Пользовательские настройки
void Set_Axis_X_Text(std::string new_text)
{
	axis_ox_text = new_text;
}

void Set_Axis_Y_Text(std::string new_text)
{
	axis_oy_text = new_text;
}

void Set_Axis_Z_Text(std::string new_text)
{
	axis_oz_text = new_text;
}

void Set_Axis_Text_Size(double new_size)
{
	axis_text_size = new_size;
}

void Set_Marks_Text_Size(double new_size)
{

	marks_text_size = new_size;
}

void Show_Z_Marks(bool boolean)
{

	build_z_marks = boolean;
}

void Set_Legend_Text(std::string* new_text)
{
	object_name = new_text;
}

void Set_Legend_Start_Pos_X(double posit)
{

	zero_x = posit;
}

void Show_Legend(bool flag)
{

	legend_need = flag;
}

void Set_Net(bool flag)
{
	net_need = flag;
}

void Set_Step_X(double new_step)
{
	step_x = new_step;
}

void Set_Step_Y(double new_step)
{
	step_y = new_step;
}

void Set_Axis_Line_Size(double new_size)
{
	axis_line_size = new_size;
}

void Set_Marks_Line_Size(double new_size)
{
	axis_mark_line_size = new_size;
}

void Set_Net_Line_Size(double new_size)
{
	net_line_size = new_size;
}

void Set_Graphic_Line_Size(double new_size)
{
	graphic_line_size = new_size;
}

void Set_Surface_Line_Size(double new_size)
{
	surface_line_size = new_size;
}

void Set_Axis_Point_Size(double new_size)
{
	axis_point_size = new_size;
}

void Set_Graphic_Point_Size(double new_size)
{
	graphic_point_size = new_size;
}

void Set_Surface_Point_Size(double new_size)
{
	surface_point_size = new_size;
}

void Set_Surface_Width(double new_size)
{
	surface_triangle_volume = new_size;
}

void Set_Text_Line_Size(double new_size)
{
	text_size = new_size;
}

void tenth_mode(bool flag)
{
	mm = flag;
}