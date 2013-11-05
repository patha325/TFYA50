#include <iostream>

using namespace std;

void hideNSeek();

int main(){

	cout << "Hello, Allan!" << endl;


    hideNSeek();
	system("pause");
	return 0;
}

void hideNSeek(){


	for (int i = 1; i<50; i++){	
		cout << i << endl;
	}
	cout << "Ready or not, here I come!" << endl;
}