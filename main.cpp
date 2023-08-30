#include <SFML/Graphics.hpp>
#include <iostream>
#include <math.h>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <exception>

using namespace sf;
using namespace std;
using Eigen::Matrix;
using MatrixXY = Matrix<float, 2, 1>; using MatrixXYZ = Matrix<float, 3, 1>;
using MatrixDistance = Matrix<float, 4, 1>;
using MatrixHeight = Matrix<float, 4, 1>;
using MatrixH = Matrix<float, 4, 2>;
using MatrixHpl = Matrix<float, 4, 3>;
using MatrixEx = Matrix<float, 2, 50>;
using MatrixExpl = Matrix<float, 3, 50>;
MatrixEx bank; MatrixExpl bankpl;
MatrixXY shot; MatrixXYZ shotpl;
Uint8 stopper;
char command;
static sf::RenderWindow window(sf::VideoMode(600, 600), "MNK Graph");

const int W = 602;
const int H = 602;
float inputX; float inputY; float inputZ;

std::random_device rd{};
std::mt19937 gen{ rd() };
std::normal_distribution<> d{ 0, 0.09 };

//Метод А - проcто добавь тригинометрии

MatrixDistance flatter(MatrixDistance dist, MatrixHeight high) {
    MatrixDistance res;
    for (int i = 0; i < 4; i++) {
        float cosin = dist(i, 0) / high(i, 0);
        float sin = sqrt(1 - pow(cosin, 2));
        res(i, 0) = dist(i, 0) * sin;
    }
    return res;
}

MatrixDistance fx(MatrixXY& x_pred, Matrix <float, 4, 2>& Sat) {
    /*
        Вектор функциональной связи между измеряемой дальностью и координатами объекта
    : param x_pred :
    : param Sat :
    : return :
     */
    MatrixDistance f;
    for (size_t i = 0; i < 4; i++) { //Количество маяков
        f(i, 0) = (sqrt(pow(Sat(i, 0) - x_pred(0, 0), 2) + pow((Sat(i, 1) - x_pred(1, 0)), 2)));
    }
    return f;
}

MatrixDistance fxpl(MatrixXYZ& x_pred, Matrix <float, 4, 3>& Sat) {
    /*
        Вектор функциональной связи между измеряемой дальностью и координатами объекта
    : param x_pred :
    : param Sat :
    : return :
     */
    MatrixDistance f;
    for (size_t i = 0; i < 4; i++) { //Количество маяков
        f(i, 0) = sqrt(pow(Sat(i, 0) - x_pred(0, 0), 2) + pow(Sat(i, 1) - x_pred(1, 0), 2) + pow(Sat(i, 2) - x_pred(2, 0), 2));
    }
    return f;
}

/*MatrixDistance fx(MatrixXY& x_pred, Matrix <float, 4, 2>& Sat) {
    /*
        Вектор функциональной связи между измеряемой дальностью и координатами объекта
    : param x_pred :
    : param Sat :
    : return :

    MatrixDistance f;
for (size_t i = 0; i < 4; i++) { //Количество маяков
    f(i, 0) = sqrt(pow(Sat(i, 0) - x_pred(0, 0), 2) + pow((Sat(i, 1) - x_pred(1, 0)), 2));
}
return f;
}*/



MatrixH Hx(MatrixXY& x_pred, Matrix <float, 4, 2>& Sat) {
    /*
        Градиентная матрица
    : param x_pred :
    : param Sat :
    : return :
       */
    Matrix <float, 4, 2> h;
    for (size_t i = 0; i < 4; i++) { //где n таковое количество маяков, сами посчитаете ты сказал
        for (size_t j = 0; j < 2; j++) { // цифра 2 потому что три координаты х, y
            h(i, j) = -(Sat(i, j) - x_pred(j, 0)) / sqrt(pow((Sat(i, 0) - x_pred(0, 0)), 2) + pow((Sat(i, 1) - x_pred(1, 0)), 2));
        }
    }
    return h;
}

MatrixHpl Hxpl(MatrixXYZ& x_pred, Matrix <float, 4, 3>& Sat) {
    /*
        Градиентная матрица
    : param x_pred :
    : param Sat :
    : return :
       */
    Matrix <float, 4, 3> h;
    for (size_t i = 0; i < 4; i++) { //где n таковое количество маяков, сами посчитаете ты сказал
        for (size_t j = 0; j < 3; j++) { // цифра 3 потому что три координаты х, y, z
            h(i, j) = -(Sat(i, j) - x_pred(j, 0)) / sqrt(pow((Sat(i, 0) - x_pred(0, 0)), 2) + pow((Sat(i, 1) - x_pred(1, 0)), 2) + pow((Sat(i, 2) - x_pred(2, 0)), 2));
        }
    }
    return h;
}

Matrix <float, 2, 1> descent_process(MatrixDistance& R, Matrix <float, 4, 2>& x_sat, MatrixXY& x_pred) {
    /*
        Итерация градиентного спуска(процесс)
        : param R : наблюдения
        : param x_sat : координаты спутников(якорей)
        : param x_prev : оценка координат объекта на k - 1 итерации
        : return : оценка координат объекта на k итерации
       */
    MatrixH h = Hx(x_pred, x_sat);
    MatrixDistance f = fx(x_pred, x_sat);
    MatrixXY x_new = x_pred + (h.transpose() * h).inverse() * h.transpose() * (R - f);
    return x_new;
}

Matrix <float, 3, 1> descent_processpl(MatrixDistance& R, Matrix <float, 4, 3>& x_sat, MatrixXYZ& x_pred) {
    /*
        Итерация градиентного спуска(процесс)
        : param R : наблюдения
        : param x_sat : координаты спутников(якорей)
        : param x_prev : оценка координат объекта на k - 1 итерации
        : return : оценка координат объекта на k итерации
       */
    MatrixHpl h = Hxpl(x_pred, x_sat);
    MatrixDistance f = fxpl(x_pred, x_sat);
    MatrixXYZ x_new = x_pred + (h.transpose() * h).inverse() * h.transpose() * (R - f);
    return x_new;
}

pair<bool, MatrixXY> mnkPsevdo(MatrixDistance& r, Matrix <float, 4, 2>& sat) {
    // Матрица Н в размерности n * 2, где n количество маяков, а столбцы отвечают за координаты х, y
    // переменная sat имеет размерность n * 2, это координаты опорных маяков
    // x это искомый вектор состояния, имеет размерность 2 * 1
    // x_pred соответственно предыдущее значение переменной Xo
    // r соответственно получаемые дальности по каждому маяку, входные данные, имеет размерность n * 1
    int n = r.rows();
    MatrixDistance R0;
    Matrix <float, 4, 2> Sat;
    for (size_t i = 0; i < n; i++) {
        if (!isnan(r(i, 0))) {
            R0(i, 0) = r(i, 0);
            for (size_t j = 0; j < 2; j++) {
                Sat(i, j) = sat(i, j);
            }
        }
    }
    MatrixXY x_pred{ {0.00001},
                      {0.00001} };

    n = R0.rows();
    MatrixDistance R;
    for (size_t i = 0; i < n; i++) {
        R(i, 0) = R0(i, 0) + d(gen);
    }

    Uint16 kkk = 0;
    float epsilon = 0.01;
    bool stop = false;
    while (true) {
        try {
            MatrixXY x_new = descent_process(R, Sat, x_pred);
            float Euc = (x_new - x_pred).norm();
            x_pred = x_new;
            bank(0, kkk) = x_pred(0, 0); bank(1, kkk) = x_pred(1, 0);
            kkk += 1;
            cout << x_pred(0, 0) << " " << x_pred(1, 0) << endl;
            //shape.setPosition(50 + x_pred(0, 0)*80, 550 - x_pred(1, 0)*80);
            if (Euc <= epsilon) {
                shot(0, 0) = x_pred(0, 0); shot(1, 0) = x_pred(1, 0);
                stopper = kkk;
                cout << "\nStoped computation" << endl;
                stop = true;
                break;
            }
            if (kkk > 50) {
                stop = false;
                break;
            }
        }
        catch (...) {
            cout << "error" << endl;
            stop = false;
            break;
        }
    }
    if (stop) {
        cout << "ok" << endl;
        return make_pair(true, x_pred);
    }
    else
        return make_pair(false, MatrixXY{});

}

pair<bool, MatrixXYZ> mnkPsevdopl(MatrixDistance& r, Matrix <float, 4, 3>& sat) {
    // Матрица Н в размерности n * 3, где n количество маяков, а столбцы отвечают за координаты х, y
    // переменная sat имеет размерность n * 2, это координаты опорных маяков
    // x это искомый вектор состояния, имеет размерность 2 * 1
    // x_pred соответственно предыдущее значение переменной Xo
    // r соответственно получаемые дальности по каждому маяку, входные данные, имеет размерность n * 1
    int n = r.rows();
    MatrixDistance R0;
    Matrix <float, 4, 3> Sat;
    for (size_t i = 0; i < n; i++) {
        if (!isnan(r(i, 0))) {
            R0(i, 0) = r(i, 0);
            for (size_t j = 0; j < 3; j++) {
                Sat(i, j) = sat(i, j);
            }
        }
    }
    MatrixXYZ x_pred{ {0.00001}, {0.00001}, {0.00001} };

    n = R0.rows();
    MatrixDistance R;
    for (size_t i = 0; i < n; i++) {
        R(i, 0) = R0(i, 0) + d(gen);
    }

    Uint16 kkk = 0;
    float epsilon = 0.0008;
    bool stop = false;
    while (true) {
        try {
            MatrixXYZ x_new = descent_processpl(R, Sat, x_pred);
            float Euc = (x_pred - x_new).norm();
            cout << Euc << endl;
            x_pred = x_new;
            bankpl(0, kkk) = x_pred(0, 0); bankpl(1, kkk) = x_pred(1, 0); bankpl(2, kkk) = x_pred(2, 0);
            kkk += 1;
            cout << x_pred(0, 0) << " " << x_pred(1, 0) << " " << x_pred(2, 0) << endl;
            if (Euc <= epsilon) {
                shotpl(0, 0) = x_pred(0, 0); shotpl(1, 0) = x_pred(1, 0); shotpl(2, 0) = x_pred(2, 0);
                stopper = kkk;
                cout << "\nStoped computation" << endl;
                stop = true;
                break;
            }
            if (kkk > 48) {
                stop = false;
                break;
            }
        }
        catch (...) {
            cout << "error" << endl;
            stop = false;
            break;
        }
    }
    if (stop) {
        cout << "ok" << endl;
        return make_pair(true, x_pred);
    }
    else
        return make_pair(false, MatrixXYZ{});

}

void Drawer(pair<bool, MatrixXY> res) {
    sf::CircleShape point(5.0); sf::CircleShape target(5.0); sf::CircleShape guess(5.0);
    point.setFillColor(sf::Color::Green); target.setFillColor(sf::Color::Red); guess.setFillColor(sf::Color::Magenta);
    point.setPosition(50, 550); window.draw(point);
    for (int i = 0; i < stopper; i++) {
        point.setPosition(50 + bank(0, i) * 80, 550 - bank(1, i) * 80);
        window.draw(point);
    }
    guess.setPosition(50.0 + res.second(0, 0) * 80, 550.0 - res.second(1, 0) * 80);
    window.draw(guess);
    target.setPosition(50.0 + inputX * 80, 550.0 - inputY * 80);
    window.draw(target);
}

void Drawerpl(pair<bool, MatrixXYZ> res) {
    sf::CircleShape point(5.0); sf::CircleShape target(5.0); sf::CircleShape guess(5.0);
    point.setFillColor(sf::Color::Green); target.setFillColor(sf::Color::Red); guess.setFillColor(sf::Color::Magenta);
    point.setPosition(50, 550); window.draw(point);
    for (int i = 0; i < stopper; i++) {
        point.setPosition(50 + bankpl(0, i) * 80, 550 - bankpl(1, i) * 80);
        window.draw(point);
    }
    guess.setPosition(50.0 + res.second(0, 0) * 80, 550.0 - res.second(1, 0) * 80);
    window.draw(guess);
    target.setPosition(50.0 + inputX * 80, 550.0 - inputY * 80);
    window.draw(target);
}


int main()
{
    bool exit = false;
    window.close();
    cout << "AMAZING PROGRAMM" << endl;
    while (!exit) {
        cout << "Please, input command (s - start, f - finish) ";
        cin >> command;
        if (command == 'f') {
            exit = true;
        }
        else if (command == 's') {
            cout << "Input X value: "; cin >> inputX;
            cout << "Input Y value: "; cin >> inputY;
            window.create(sf::VideoMode(600, 600), "MNK Graph");
            Matrix <float, 4, 2> Sat{ {1, 10.5}, {1, 2}, {5, 11.5}, {5, 2} };
            MatrixXY point1{ {inputX},{inputY} };
            MatrixDistance r = fx(point1, Sat);
            pair<bool, MatrixXY> result = mnkPsevdo(r, Sat);
            if (result.first) {
                float x = result.second(0, 0);
                float y = result.second(1, 0);
                cout << "x = " << x << ", y = " << y << endl;
            }
            else {
                cout << "Something went wrong!" << endl;
                window.close();
            }
            while (window.isOpen())
            {
                sf::Event event;
                while (window.pollEvent(event))
                {
                    if (event.type == sf::Event::Closed) {
                        window.close();
                    }
                }
                window.clear(sf::Color::White);
                Drawer(result);
                window.display();
            }
        }
        else if (command == 'p') {
            cout << "Input X value: "; cin >> inputX;
            cout << "Input Y value: "; cin >> inputY;
            cout << "Input Z value: "; cin >> inputZ;
            window.create(sf::VideoMode(600, 600), "MNK Graph");
            Matrix <float, 4, 3> Satpl{ {1, 10.5, 2.0}, {1, 2, 2.0}, {5, 11.5, 2.0}, {5, 2, 2.0} };
            MatrixXYZ point2{ {inputX},{inputY},{inputZ} };
            MatrixDistance rpl = fxpl(point2, Satpl);
            pair<bool, MatrixXYZ> result = mnkPsevdopl(rpl, Satpl);
            if (result.first) {
                float x = result.second(0, 0);
                float y = result.second(1, 0);
                float z = result.second(2, 0);
                cout << "x = " << x << ", y = " << y << ", z = " << z << endl;
            }
            else {
                cout << "Something went wrong!" << endl;
                window.close();
            }
            while (window.isOpen())
            {
                sf::Event event;
                while (window.pollEvent(event))
                {
                    if (event.type == sf::Event::Closed) {
                        window.close();
                    }
                }
                window.clear(sf::Color::White);
                Drawerpl(result);
                window.display();
            }
            //bankpl.setZero();
        }
        else {
            cout << "Sorry, wrong input" << endl;
        }
    }
    return 0;
}