using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;



namespace LibraryCorrelations
{
    class add_on_function

    {   public
        add_on_function()
        {
        }
        
        	public int electrode_find_file_position (string [] name_electrodes_file, string key)
            {

	  		for (int i=0; i < name_electrodes_file.Length; i++)
	   	 	{
	      	if (Path.GetFileNameWithoutExtension((name_electrodes_file[i]).Substring(name_electrodes_file[i].LastIndexOf("_") + 1))==key)
			return i;
	    	}
	 	 	return -1;
			}
    

    	public bool contains (int [] vec, int key)
       		{
            for (int i=0; i< vec.Length; i++)
            {
                if (vec[i] == key)
                    return true;
            }
            return false;
          
      		 }

     	public bool contains(string[] vec, string key)
       		 {
            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i] == key)
                    return true;
            }
            return false;

       		 }


        public char[] nonzeros(char[] vec)
        {
            int j = 0;


            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i] != 0)
                    j++;

            }

            char[] n_zeros = new char[j];

                j = 0;

            for (int i = 0; i < vec.Length; i++)
            {

                if (vec[i] != 0)
                {
                    n_zeros[j] = vec[i];
                    j++;
                }
            }

            return n_zeros;
        }

        public double[] nonzeros(double[] vec)
        {
            int j = 0;

            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i] != 0)
                    j++;

            }

            double[] n_zeros = new double[j];

            j = 0;

            for (int i = 0; i < vec.Length; i++)
            {

                if (vec[i] != 0)
                {
                    n_zeros[j] = vec[i];
                    j++;
                }
            }

            return n_zeros;
        }

        public int[] nonzeros(int[] vec)
        {
            int j = 0;


            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i] != 0)
                    j++;

            }

            int[] n_zeros = new int[j];

            j = 0;

            for (int i = 0; i < vec.Length; i++)
            {

                if (vec[i] != 0)
                {
                    n_zeros[j] = vec[i];
                    j++;
                }
            }

            return n_zeros;
        }

        public int find (int [] vec, int key)
        {
           for (int i = 0; i < vec.Length; i++)
           {
               if (vec[i] == key)
                   return i;
           }
           return -1;

        }

        public int find(string[] vec, string key)
        {
            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i] == key)
                    return i;
            }
            return -1;
        }

        public int find(double[] vec, double key)
        {
            for (int i = 0; i < vec.Length; i++)
            {
                if (vec[i] == key)
                    return i;
            }
            return -1;

        }


    



	
    }
}

