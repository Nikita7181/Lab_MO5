#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include "supportLib.hpp"
#include "pbPlots.hpp"

#define pi 3.1415
#define K 100
#define a 25
#define P 0.95
#define epsilon 0.01
#define r 3
#define M (r - 1) / 2

double get_omega(std::vector<double> f)
{
    int size = f.size();
    double sum = 0.0;
    for (int k = 1; k < size; k ++)
    {
        sum = sum + pow((f[k]-f[k-1]),2);
    }
    return sqrt(sum);
}

double get_delta(std::vector<double> f_filtered, std::vector<double> f_noised)
{
    int size = f_filtered.size();
    double sum = 0.0;

    for (int i = 0; i < size; i++)
    {
        sum = sum + pow(f_filtered[i] - f_noised[i], 2);
    }

    return sqrt(sum/size);
}

double get_dist (double delta, double omega)
{
    return sqrt(pow(delta, 2) + pow(omega, 2));
}

double filter (std::vector<double> f_noised, std::vector<double> alpha, int k)
{
    double sum = 0.0;
    for (int j = k - M -1; j <= k + M -1; j++)
    {
        sum = sum + pow(f_noised[j], 2) * alpha[j + M + 1 - k];
    }
    return sqrt(sum);

}

std::vector <double> generate_alpha(int window)
{
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0,  1);
    double alpha_centr =  distribution(generator);
    std::vector<double> result; result.resize(window, 0);
    result[M] = alpha_centr; // по формуле M+1, но нумерация с 0
    double sum = 0.0;
    for (int m = 2; m <= M; m++)
    {
        for (int s = m; s < r - m; s++)
        {
            sum =  sum + result[s];
        }
        std::uniform_real_distribution<double> distribution2(0,  1 - sum);
        result[m-1] = result[r-m] = 0.5 * distribution2(generator);
        sum = 0.0;
    }


    sum =  0.0; //result[M];
    for (int s = 2; s <= r - 1; s++)
        sum = sum + result[s-1];

    result[0] = result[r-1] =  0.5 * ( 1 - sum);
    sum = 0.0;
    for (int i = 0; i < result.size(); i++)
        sum = sum + result[i];
    for (int i = 0; i < result.size(); i++)
        result[i] = result[i] / sum;
    return result;
}

void print_v (std::vector<double> v, int start, int end = 0)
{
    if (end == 0)
        end=v.size()-1;
    for (int i = start; i <= end; i++)
        if (v[i] < 0 )
            std::cout << std::fixed << std::setprecision(5)  << v[i] << "  ";
        else
            std::cout << std::fixed << std::setprecision(5) << " " << v[i] << "  ";
    std::cout << "\n";
}

int main() {
    double x_min = 0;
    double x_max = pi;

    std::vector<double> x;
    std::vector<double> f;
    std::vector<double> noise;
    std::vector<double> f_noised;
    std::vector<std::vector<double>> alpha;
    std::vector<double> f_filterd;

    for (int k = 0; k <= K; k++)
    {
        double x_k = x_min + k * (x_max - x_min) / K;
        x.push_back(x_k);
        f.push_back(sin(x_k) + 0.5);

        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distr(-1 * a, a);
        double dice = distr(generator);

        noise.push_back(0.0 + dice / 100.0);
        f_noised.push_back( f[k] + noise[k]);

    }
    f_filterd.resize(f.size(),0);

    //print_v (f, 0);
    //print_v (noise, 0, K);
    //print_v (f_noised, 0, K);


    int lambda_size = 11;
    int alpha_size = r;
    alpha.resize(lambda_size, std::vector<double> (alpha_size, 0));

    int expreriment_count = static_cast<int> (floor( log(1-P)/log (1 - epsilon/(x_max - x_min)))); //N

    std::vector<std::vector<double>> result_table(lambda_size, std::vector<double> (5, 0));
    double lambda = 0;
    for (int step = 0; step < lambda_size ; step++)
    {
        std::vector<double> best_alpha; best_alpha.resize(r, 0);
        double min_omega = 9999;
        double min_delta = 9999;
        double min_criteria = 99999;
        for (int ex = 0; ex < expreriment_count; ex++)
        {
            alpha[step] = generate_alpha(r);

            for(int k = M+1 ; k < K; k++)
            {
                f_filterd[k]= filter(f_noised,alpha[step],k);
            }
            double tmp_omega = get_omega(f_filterd);
            double tmp_delta = get_delta(f_filterd, f_noised);
            double tmp_criteria = lambda * tmp_omega + (1 - lambda) * tmp_delta; //J  //??????????????get_dist(tmp_delta, tmp_omega);//
            if (tmp_criteria < min_criteria)
            {
                min_criteria = tmp_criteria;
                min_omega = tmp_omega;
                min_delta = tmp_delta;
                best_alpha = alpha[step];
            }
        }
        alpha[step] = best_alpha;


        std::cout << "lambda = " << lambda << " J = " << min_criteria << " omega = " << min_omega << " delta = " << min_delta << " alpha: [";
        for (int i = 0; i < r; i++)
            std::cout << "  " << alpha[step][i] ;
        double dist = get_dist (min_delta, min_omega);
        std::cout << " ] dist = " << dist <<"\n";
        //std::cout << " ]\n";
        result_table[step][0] = lambda;
        result_table[step][1] = min_criteria;
        result_table[step][2] = min_omega;
        result_table[step][3] = min_delta;
        result_table[step][4] = get_dist(min_delta, min_omega) ;////lambda * min_omega + (1 - lambda) * min_delta;;

        lambda = lambda + 0.1;
    }

    int index = 0;
    for (int i = 1; i < lambda_size; i++)
    {
        if (result_table[index][4] > result_table[i][4])
            index = i;
    }

    std::cout << "\n best = " << index + 1 << "\n";
    std::cout << " lambda*      J       omega     delta      dist\n";

    print_v(result_table[index],0);

    std::cout << " alpha :";
    print_v(alpha[index],0);

    for(int k = 0; k < K; k++)
    {
        f_filterd[k]= filter(f_noised,alpha[index],k);
    }

    //print_v (f_noised, 0, K);
    //print_v (f, 0, K);
    //print_v (f_filterd, 0, K);


    RGBABitmapImageReference * imageRef = CreateRGBABitmapImageReference();
    StringReference *errorMessage = new StringReference();

    ScatterPlotSeries * f_series = GetDefaultScatterPlotSeriesSettings();
    f_series->xs = &x;
    f_series->ys = &f;
    f_series->linearInterpolation = false;
    f_series->pointType = toVector(L"dots");
    f_series->lineType = toVector(L"dashed");
    f_series->lineThickness = 2;
    f_series->color = CreateRGBColor(1, 0, 0);;

    ScatterPlotSeries * f_noise_series = GetDefaultScatterPlotSeriesSettings();
    f_noise_series->xs = &x;
    f_noise_series->ys = &f_noised;
    f_noise_series->linearInterpolation = true;
    f_noise_series->pointType = toVector(L"dots");
    f_noise_series->lineType = toVector(L"dashed");
    f_noise_series->lineThickness = 2;
    f_noise_series->color = CreateRGBColor(0, 1, 0);

//    ScatterPlotSeries * noise_series = GetDefaultScatterPlotSeriesSettings();
//    noise_series->xs = &x;
//    noise_series->ys = &noise;
//    noise_series->linearInterpolation = true;
//    noise_series->pointType = toVector(L"dots");
//    noise_series->lineType = toVector(L"dashed");
//    noise_series->lineThickness = 2;
//    noise_series->color = CreateRGBColor(0, 1, 1);

    ScatterPlotSeries * f_filtered_series = GetDefaultScatterPlotSeriesSettings();
    f_filtered_series->xs = &x;
    f_filtered_series->ys = &f_filterd;
    f_filtered_series->linearInterpolation = true;
    f_filtered_series->pointType = toVector(L"dots");
    f_filtered_series->lineType = toVector(L"dashed");
    f_filtered_series->lineThickness = 2;
    f_filtered_series->color = CreateRGBColor(0, 0 , 1);

    ScatterPlotSettings * settings = GetDefaultScatterPlotSettings();
    settings->width = 2000;
    settings->height = 700;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    if(r == 3)
        settings->title = toVector(L"Plots of the original signal, noise, cleaned signal for r = 3");
    if(r == 5)
        settings->title = toVector(L"Plots of the original signal, noise, cleaned signal for r = 5");
    settings->xLabel = toVector(L"X axis");
    settings->yLabel = toVector(L"Y axis");
    settings->scatterPlotSeries->push_back(f_series);
    settings->scatterPlotSeries->push_back(f_noise_series);
//    settings->scatterPlotSeries->push_back(noise_series);
    settings->scatterPlotSeries->push_back(f_filtered_series);

    DrawScatterPlotFromSettings(imageRef, settings, errorMessage);

    std::vector<double> * pngData = ConvertToPNG(imageRef->image);
    WriteToFile(pngData, "../Output.png");
    DeleteImage(imageRef->image);

    RGBABitmapImageReference * imageRef2 = CreateRGBABitmapImageReference();

    ScatterPlotSettings *settings2 = GetDefaultScatterPlotSettings();
    settings2->width = 1000;
    settings2->height = 700;
    settings2->autoPadding = true;
    settings2->autoBoundaries = false;
    settings2->xMin = 0.0;
    settings2->xMax = 2.5;
    settings2->yMin = 0.0;
    settings2->yMax = 0.2;
    settings2->showGrid = true;

    settings2->title = toVector(L"Graphic display of the found approximations to the optimal criteria");
    settings2->xLabel = toVector(L"omega axis");
    settings2->yLabel = toVector(L"delta axis");

    std::vector<std::vector<double>> w_d_ys1;
    w_d_ys1.resize(result_table.size() + 1, std::vector<double>(1,0));
    std::vector<std::vector<double>> w_d_xs1;
    w_d_xs1.resize(result_table.size() + 1, std::vector<double>(1,0));

    for(int i = 0; i < result_table.size(); i++)
    {
        w_d_xs1[i][0] = result_table[i][2];
        w_d_ys1[i][0] = result_table[i][3];
    }

    w_d_ys1[11][0] = 0.0000000001;
    w_d_xs1[11][0] = 0.0000000001;

    std::vector<ScatterPlotSeries*> w_d;
    w_d.resize(result_table.size()+1);
    for(int i = 0; i < result_table.size() + 1; i++)
    {
        w_d[i] = new ScatterPlotSeries();

        w_d[i]->linearInterpolation  = false;
        w_d[i]->lineType = toVector(L"solid");
        w_d[i]->lineThickness = 1.0;
        w_d[i]->xs = &w_d_xs1[i];
        w_d[i]->ys = &w_d_ys1[i];
        w_d[i]->color = GetBlack();
        w_d[i]->color = CreateRGBColor((random()%100)/ 100.0,  (random()%70)/ 100.0,  (random()%50)/ 100.0);
    }

    w_d[0]->pointType= toVector(L"triangles");
    w_d[1]->pointType= toVector(L"dots");
    w_d[2]->pointType= toVector(L"crosses");
    w_d[3]->pointType= toVector(L"circles");
    w_d[4]->pointType= toVector(L"filled triangles");
    w_d[5]->pointType= toVector(L"triangles");
    w_d[6]->pointType= toVector(L"dots");
    w_d[7]->pointType= toVector(L"crosses");
    w_d[8]->pointType= toVector(L"circles");
    w_d[9]->pointType= toVector(L"filled triangles");
    w_d[10]->pointType =toVector(L"triangles");
    w_d[11]->pointType =toVector(L"dots");
    for(int i = 0; i < result_table.size() + 1; i++) {
        settings2->scatterPlotSeries->push_back(w_d[i]);
    }

    DrawScatterPlotFromSettings(imageRef2, settings2, errorMessage);

    std::vector<double> * pngData2 = ConvertToPNG(imageRef2->image);
    WriteToFile(pngData2, "../Output2.png");
    DeleteImage(imageRef2->image);
}
