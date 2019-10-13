#include <bits/stdc++.h>
using namespace std;

int main(){
	fstream file;
	string word, filename;

//	filename = "./TRANSMIT_Be_Window_100.txt";
	filename = "./TRANSMIT_Be_Window_10000.txt";

	file.open(filename.c_str());

	int dataflag = 0;
	int datacol = 0;

	while (file >> word){
//		cout << word << endl;
		if (word == "T") dataflag = 1;
		if (dataflag == 1){
			++datacol;
			if (datacol%10 == 1){
				string head = word;
			}else if (datacol%10 == 2){
				int Nion = stoi(word);
			}else if (datacol%10 == 3){
				int A = stoi(word);
			}else if (datacol%10 == 4){
				double energy = stod(word);
			}else if (datacol%10 == 5){
				double depth = stod(word);
			}else if (datacol%10 == 6){
				double LposY = stod(word);
			}else if (datacol%10 == 7){
				double LposZ = stod(word);
			}else if (datacol%10 == 8){
				double cosx = stod(word);
			}else if (datacol%10 == 9){
				double cosy = stod(word);
			}else{
				double cosz = stod(word);
			}
		}
	}

	cout << datacol/10 << "ions analyzed" << endl;

	return 0;
}
