#include <iostream>

using namespace std;

void hideNSeek();

int main(){

	cout << "Hello, World!" << endl;


    hideNSeek();
	system("pause");
	return 0;
}

void hideNSeek(){

	for (int i = 0; i<10; i++){
	
		cout << i << endl;
	}
	cout << "Ready or not, here I come!" << endl;
}