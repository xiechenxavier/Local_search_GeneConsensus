#include <iostream>
#include <algorithm>
#include <vector>
 
int main()
{
    using namespace std;
 
    vector<int> vec;
 
    vec.push_back(0);
    vec.push_back(1);
    vec.push_back(2);
    vec.push_back(3);
    vec.push_back(4);
    vec.push_back(5);
    vec.push_back(6);
 
    vector<int>::iterator it = find(vec.begin(), vec.end(), 6);
    // vec.insert(vec.begin() + 2, vec[4]);
    // vec.erase(vec.begin() + 5); //1,2,5,3,4,6
 
    vec.insert(vec.begin() + 5 + 1, vec[1]);
	vec.erase(vec.begin() + 1);

    for(auto i=0;i<vec.size();i++){
        cout<<vec[i]<<endl;
    }

    if (it != vec.end())
        cout<<*it<<endl;
    else
        cout<<"can not find"<<endl;
 
    return 0;
}