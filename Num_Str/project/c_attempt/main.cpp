#include <array>

const int N=120;
const int t_end=20000;



class system
{
	double step = 0;
	double dx = 1/(N-1);
	double dt = 1e-4;
	double r_e = 0.999;

	std::array<std::array<std::array<double, N>, 3>, t_end> save;
	std::array<std::array<double, N>, 3> Abl_tm1;

	void derrivative();
	void AB2();
};

void system::AB2()
{
	save[step] = Abl_tm1;
	derrivative();
	save[step] = Abl_tm1;
	step++;
	return;
}

int main(void)
{

	return 0;
}
