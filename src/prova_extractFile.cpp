#include <iostream>
#include <iomanip>
#include <fstream>

int main()
{
    // Number of data
    int n = 256; 
    std::fstream myfile("Q.txt", std::ios_base::in);
    float Qarray [n][n];
    for (int i = 0; i<n; i++ )
        for (int j = 0; j<n; j++ )
        {
            {
                std::string a;
                myfile >> a;
                // change string to float
                float b = std::stof (a);
                Qarray[i][j] = b;
            }
        }
    std::cout << Qarray[255][255] << std::endl;        
    return 0;  
}
