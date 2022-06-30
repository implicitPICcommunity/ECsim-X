
#ifndef ALLOC_H
#define ALLOC_H

#include <stdlib.h>
#include <stdio.h>

/*! The allocator for 4D array */
template < class type > type **** _new_4_array(int sz1, int sz2, int sz3, int sz4) {

  type ****all_x;
  type ***all_y;
  type **all_z;
  type *all_r;

  all_x = new type ***[sz1];
  all_y = new type **[sz1 * sz2];
  all_z = new type *[sz1 * sz2 * sz3];
  all_r = new type[sz1 * sz2 * sz3 * sz4];

  type ****result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
      for (int k = 0; k < sz3; k++, all_r += sz4) {
        result[i][j][k] = all_r;
      }
    }
  }

  return result;
}

/*! The assigment for 4D array */
template < class type > type **** _assign_4_array(int sz1, int sz2, int sz3, int sz4, type **** org) {

  type ****all_x;
  type ***all_y;
  type **all_z;
  type *all_r;

  all_x = org;
  all_y = org[0];
  all_z = org[0][0];
  all_r = org[0][0][0];

  type ****result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
      for (int k = 0; k < sz3; k++, all_r += sz4) {
        result[i][j][k] = all_r;
      }
    }
  }

  return result;
}

/*! Deallocator for 4D arrays */
template < class type > void delArr4(type **** arr, int dummyx, int dummyy, int dummyz) {
  delete[]arr[0][0][0];
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 3D array */
template < class type > type *** _new_3_array(int sz1, int sz2, int sz3) {

  type ***all_x;
  type **all_y;
  type *all_z;

  all_x = new type **[sz1];
  all_y = new type *[sz1 * sz2];
  all_z = new type[sz1 * sz2 * sz3];

  type ***result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
      if (i == 0 and j == 0){
	 //printf("i = 0 e j = 0 all_z \n");     
         //printf("_new_3_array, all_z = %p \n", all_z);
      }
    }
  }

  //uncommet this to verify the mapping from 3D to 1D array
  //printf("_new_3_array, result[0][0] = %p \n", result[0][0]);
  //printf("_new_3_array, result = %p \n", result);
  return result;

}

/*! The assignment for 3D array */
template < class type > type *** _assign_3_array(int sz1, int sz2, int sz3, type *** org) {

  type ***all_x;
  type **all_y;
  type *all_z;

  all_x = org;
  all_y = org[0];
  all_z = org[0][0];

  type ***result = org;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
    }
  }
  
  return result;
}

/*! Deallocator for 3D arrays */
template < class type > void delArr3(type *** arr, int dummyx, int dummyy) {
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 2D array */
template < class type > type ** _new_2_array(int sz1, int sz2) {

  type **all_x;
  type *all_y;

  all_x = new type *[sz1];
  all_y = new type[sz1 * sz2];

  type **result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
  }

  return result;

}

/*! The assignment for 2D array */
template < class type > type ** _assign_2_array(int sz1, int sz2, type ** org) {

  type **all_x;
  type *all_y;

  all_x = org;
  all_y = org[0];

  type **result = org;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
  }

  return result;
}

/*! Deallocator for 2D arrays */
template < class type > void delArr2(type ** arr, int dummyx) {
  delete[]arr[0];
  delete[]arr;
}

struct A5D {

  A5D(int n1, int n2, int n3, int n4, int n5) :
    l1(n2*n3*n4*n5), l2(n3*n4*n5), l3(n4*n5),l4(n5),
    r2(n1), r3(n1*n2), r4(n1*n2*n3), r5(n1*n2*n3*n4)
  {
    data = (double*) calloc(n1*n2*n3*n4*n5,sizeof(double));
    cudaMemPrefetchAsync(data, sizeof(double)*n1*n2*n3*n4*n5, 0, 0);
  }

  double& l(int i1, int i2, int i3, int i4, int i5) {
    return data[l1*i1+i2*l2+i3*l3+i4*l4+i5];
  }
  double& r(int i1, int i2, int i3, int i4, int i5) {
    return data[i1+i2*r2+i3*r3+i4*r4+i5*r5];
  }
 
  void delArr5(){
       free(data);
  }

  const int l1, l2, l3, l4;
  const int r2, r3, r4, r5;

  double *data;

};


#define newArr4(type,sz1,sz2,sz3,sz4) _new_4_array<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3(type,sz1,sz2,sz3) _new_3_array<type>((sz1),(sz2),(sz3))
#define newArr2(type,sz1,sz2) _new_2_array<type>((sz1),(sz2))

#define asgArr2(type,sz1,sz2,org) _assign_2_array<type>((sz1),(sz2),(org))
#define asgArr3(type,sz1,sz2,sz3,org) _assign_3_array<type>((sz1),(sz2),(sz3),(org))
#define asgArr4(type,sz1,sz2,sz3,sz4,org) _assign_4_array<type>((sz1),(sz2),(sz3),(sz4),(org))

#endif
