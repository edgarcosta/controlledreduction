// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

//saves in a very raw way everything into a text file
bool de_Rham_local::save(const char * filename)
{
    if(verbose)
        cout <<"Saving de_Rham_local to: "<< filename << endl;

    ofstream fileswp;
    char fileswpname[sizeof(filename)+4];
    strcpy(fileswpname, filename);
    strcat(fileswpname,".swp");
    fileswp.open(fileswpname);

    int64_t i;
    ofstream file;
    Vec<zz_p> fbar_vector;
    map< Vec<int64_t>, zz_p, vi64less>::const_iterator it;
    pair< Vec<int64_t> , Mat<ZZ_p> > *solve_pair;

    fbar_vector.SetLength( tuple_list[d].length() );
    for( it = fbar.begin(); it != fbar.end(); it++)
        fbar_vector[ tuple_dict[d][it->first] ] = it->second;

    file.open(filename, ios::out | ios::trunc);

    

    if(file.is_open() && fileswp.is_open())
    {
        file << p <<endl;
        file << precision <<endl;
        file << verbose <<endl;
        file << save_memory << endl;
        file << n <<endl;
        file << d <<endl;
        file << Hilbert_J <<endl;
        file << fbar_vector << endl;

        for( i = 0; i < n ; i++)
        {
            solve_pair = get_solve_J( (i+1)*d - (n+1) );
            file << solve_pair->first <<endl;
            file << solve_pair->second <<endl;
        }
        solve_pair = get_solve_J((d - 2)*(n + 1) + 1);
        file << solve_pair->first <<endl;
        file << solve_pair->second <<endl;
        file.close();
        fileswp.close();
        remove(fileswpname);
        return true;
    }
    else
    {
        cout << "Wasn't able to open "<<filename<<" or "<<fileswpname<<endl;
        return false;
    }
}
