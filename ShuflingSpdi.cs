using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.IO;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using LibraryCorrelations;
using CenterSpace.NMath.Core; //esto es para hacer suavizados

namespace  ShuflingSpdi
{
  	class shufling_spdi
  	{
    
		public shufling_spdi(){}

        
        public void create_dir_dt(string name_experiment)
    	{     
                string path = name_experiment;

                try
                {
                    // Determine whether the directory exists.
                    if (Directory.Exists(path))                    
                        Directory.Delete(path,true); 
                        //Console.WriteLine("Eliminando: {0}",path);
                    // Try to create the directory.
                    DirectoryInfo di = Directory.CreateDirectory(path);
                    Console.WriteLine(path);
                }

                catch (Exception e)
                {
                    Console.WriteLine("El proceso falló: {0}", e.ToString());
                }
                finally {}
                
               
		}



        public void limpiar(string[] address_parent,int numero_de_experimentos,string path_delay,string path_inhib)
        {              
            string address = "o";
			
			int number_of_experimental_set = address_parent.Length;
		
            for (int num_set = 0; num_set < numero_de_experimentos; num_set++)
            {                               
                address = address_parent[num_set].ToString();
 
                create_dir_surr(address); 
            }

            string[] add = Directory.GetFiles(path_delay,"*");
                       
            if(add.Length != 0)
            {   
                 

                string[]  FileName = Directory.GetFiles(path_delay,"*");

                for(int ix =0; ix < FileName.Length; ix++)
                {
                    if (File.Exists(FileName[ix]))
                    {
                    File.Delete(FileName[ix]);
                    }
                    
                }
            }   

            string[] xadd = Directory.GetFiles(path_inhib,"*");
                       
            if(xadd.Length != 0)
            {   
               

                string[]  FileName = Directory.GetFiles(path_inhib,"*");

                for(int ix =0; ix < FileName.Length; ix++)
                {
                    if (File.Exists(FileName[ix]))
                    {
                    File.Delete(FileName[ix]);
                    }
                    
                }
            }  

        }

        void create_dir_surr(string name_experimento)
    	{
            string[] addres = Directory.GetDirectories(name_experimento,"*surr");
            
            if (addres.Length !=0)
            {  
                string path = addres[0].ToString();

                try
                {   
                    // Determine whether the directory exists.
                    if (Directory.Exists(path))    
                        Directory.Delete(path,true); 
                        //Console.WriteLine("Eliminando: {0}",path);
                    // Try to create the directory.
                    //DirectoryInfo di = Directory.CreateDirectory(path);

                }

                catch (Exception e)
                {
                    Console.WriteLine("El proceso falló: {0}", e.ToString());
                }
                finally {}
                
            }   
            
            string[] adx = Directory.GetDirectories(name_experimento,"*sua*");
            string[] zdx = Directory.GetDirectories(adx[0],"ncc*");

            if (zdx.Length !=0)
            {  
                string pxath = zdx[0].ToString();
                 
                try
                {   
                    
                    if (Directory.Exists(pxath))    
                        Directory.Delete(pxath,true); 
       
                  

                }

                catch (Exception e)
                {
                    Console.WriteLine("El proceso falló: {0}", e.ToString());
                }
                finally {}
                
            }
            
            

		}
         
       	void create_dir_files(string name_experimento)
    	{
            string[] addres = Directory.GetDirectories(name_experimento,"ncc*");
            
            if (addres.Length !=0)
            {  
                string path = addres[0].ToString();

                try
                {
                    // Determine whether the directory exists.
                    if (Directory.Exists(path))                    
                        Directory.Delete(path,true); 

                    // Try to create the directory.
                    //DirectoryInfo di = Directory.CreateDirectory(path);

                }

                catch (Exception e)
                {
                    Console.WriteLine("El proceso falló: {0}", e.ToString());
                }
                finally {}
                
            }   
		}
        void filesDelete(string address)
        {            
            string[] add = Directory.GetFiles(address,"*.bin");
                       
            if(add.Length != 0)
            {   
                 

                string[]  FileName = Directory.GetFiles(address,"*.bin");

                for(int ix =0; ix < FileName.Length; ix++)
                {
                    if (File.Exists(FileName[ix]))
                    {
                    File.Delete(FileName[ix]);
                    }
                    
                }
            }   
             
        } 

        public void shufling(bool fsuavz,string[] address_parent,int numero_de_experimentos, float fs, float mfr, float bin, float window, int num_surr, int wsurrgte,string path_mcc,string path_delay,string path_inhib, bool delay_inhib,bool recompute_inhib)
        {              
            string address = "o";
            int[] id_el={0};
			int isbin=0;	
            double n_thresh =0.0d;
            
			//donde, fs: creo que es la frecuencia de muestreo, window: es el ancho de la ventana, bin: es el ancho de bin
			// mfr = spikes/s=0.1:frecuencia de disparo.
					
			int number_of_experimental_set = address_parent.Length;
		
            for (int num_set = 0; num_set < numero_de_experimentos; num_set++)
            {                               
                address = address_parent[num_set].ToString();

                Console.WriteLine(address);
                
                create_dir_surr(address);  

                string[] xmadres = Directory.GetDirectories(address,"Exp*");

                address = xmadres[0];
               
                create_dir_files(address);
            
                string name = Path.GetFileName(address);

                string[] namefile2 = Directory.GetFiles(address);

                string[] namefile = new string[namefile2.Length];


                int number_electrodes = namefile.Length;

                // It's convenient to sort the name of the electrodes' file in crescent order, in order to semplify
                //the computing of the Connectivity Matrix (CM) and to better order the cross and partia correlogram

                int[] id_electrodes = new int[namefile.Length];

                // At first we save all the electrodes indexes in the order that they are red by the software

                for (int i = 0; i < number_electrodes; i++){
                    id_electrodes[i] = int.Parse(Path.GetFileNameWithoutExtension((namefile2[i]).Substring(namefile2[i].LastIndexOf("_") + 1)));
                    }
                // Then we sort the array in crescent order, by using a built-in function

                Array.Sort(id_electrodes);

                // Now we use a function that simply returns the position of the electrodes file with a specific index (id)

                id_el = id_electrodes;

                add_on_function to_sort_electrodes = new add_on_function();

                for (int i = 0; i < namefile.Length; i++)

                //In this way we built a string vector with the electrodes sorted in crescent order, and we will use this vector in the rest
                //of the code
                {
                    int file_position = to_sort_electrodes.electrode_find_file_position(namefile2, id_electrodes[i].ToString());

                    if (file_position == -1)

                        namefile[i] = namefile2[i];

                    else

                        namefile[i] = namefile2[file_position];
                }

                FileStream electrode = new FileStream(namefile[0], FileMode.Open);
                int n_samples = 0;

                if (namefile[0].Contains(".bin"))
                {
                    isbin = 1;
                    BinaryReader El = new BinaryReader(electrode);
                    n_samples = El.ReadInt32();
                    El.Dispose();
                }

                else

                {
                    isbin = 0;
                    StreamReader El = new StreamReader(electrode);
                    n_samples = int.Parse(El.ReadLine());
                    El.Dispose();
                }


                int n_files = namefile.Length;
                electrode.Close();

                int temp = 0;
                int dim_total_spikes = 0;

                //array containing the index of last spike for every electrode in the cumulative array of spike times
                int[] dim_electrode2 = new int[number_electrodes];

                //#endregion


                //We compute the correlation window's size

                float L = (float)Math.Truncate((window / bin + 0.5));


                // #region reading files


                LibraryCorrelations.add_on_function funzione = new LibraryCorrelations.add_on_function();
                int[] non_silent = new int[number_electrodes];
                int[] silenti = new int[number_electrodes];

                //array containing the number of spikes of each electrode
                int[] single_dimension = new int[number_electrodes];

                for (int k = 0; k < number_electrodes; k++)
                {
                    FileStream electrode_k = new FileStream(namefile[k], FileMode.Open);

                    // check if the file is binary of text

                    if (isbin == 0)

                    {

                        StreamReader El_k = new StreamReader(electrode_k);

                        // we determine the number of spikes (i.e. the number of elements contained in each file) 

                        El_k.ReadLine();
                        //the second element in every file is the number of spikes
                        temp = int.Parse(El_k.ReadLine());
                        El_k.Dispose();

                    }

                    else

                    {

                        BinaryReader El_k = new BinaryReader(electrode_k);

                        // we determine the number of spikes (i.e. the number of elements contained in each file) 

                        El_k.ReadInt32();

                        //the second element in every file is the number of spikes
                        temp = El_k.ReadInt32();
                        El_k.Dispose();
                    }
                    //build non-silent and silent electrode indexes and build vector of single and cumulative dimensions

                    //build non-silent and silent electrode indexes and build vector of single and cumulative dimensions

                    double temp_fr = 0;
                    if (temp / (n_samples / fs) >= mfr)

                    {
                        non_silent[k] = k + 1;
                        single_dimension[k] = temp;
                        dim_electrode2[k] = dim_total_spikes + temp;
                        //update the cumulative dimension
                        dim_total_spikes = dim_total_spikes + temp;

                        //new 10/09/2016 -> compute mfr
                        temp_fr = (double)(temp / (n_samples / fs));
                
                    }

                    else
                        silenti[k] = k + 1;


                    electrode_k.Close();
                }

                Console.WriteLine("numero de muestas : {0}",n_samples);
                dim_electrode2 = funzione.nonzeros(dim_electrode2);
                single_dimension = funzione.nonzeros(single_dimension);


                //array containg the (occupied)bin times for the whole set of electrodes
                //i.e. if the bin of value -n indicates that a spike occured in the interval [n*(1/fs), (n+1)*(1/fs)]

                int[] spikes_data_2 = new int[dim_total_spikes];

                #region Channel Division

                //We use the minimum mfr set by the intial graphic interface to determine the more active 
                //channels (i.e. the non-silent channels); only these electrodes are involved in the correlation 
                //analysis, the other ones (i.e. the silent channels) are discarded. 
                non_silent = funzione.nonzeros(non_silent);
                silenti = funzione.nonzeros(silenti);

                for (int i = 0; i < non_silent.Length; i++)
                    --non_silent[i];
                for (int i = 0; i < silenti.Length; i++)
                    --silenti[i];


                int nchNS = non_silent.Length;

                int[] dim_electrode = new int[nchNS];

                //array containg the index of the first spike for every electrode in the cumulative array of spikes
                int[] dim_electrode3 = new int[nchNS];

                for (int i = 1; i < nchNS; i++){
                    dim_electrode3[i] = dim_electrode2[i - 1];}

                    #endregion


                    #region surrogate
                    
                    string[] addr = new string[num_surr + 1];
                    addr[0] = address;

                    #region if the files are txt
                if (isbin == 0)

                {
                    var seed = Environment.TickCount;
                    int w = (int)(wsurrgte * (fs / 1000));
                    //addr[0] = address;
                    var rand = new Random(seed);
                    int dim;
                    FileStream[] eee = new FileStream[num_surr];
                    StreamWriter[] sss = new StreamWriter[num_surr];
                    string first_father = Directory.GetParent(namefile[0]).ToString();
                    //Console.WriteLine(first_father);
                    string name_phase = first_father.Substring(first_father.LastIndexOf("/"));
                    
                    
                    //first_father = Directory.GetParent(first_father).ToString();
                    //first_father = Directory.GetParent(first_father).ToString();

                    ///fix////
                    string name_exp = first_father.Substring(first_father.LastIndexOf("/"));
                    
                    string address_father = string.Concat(Directory.GetParent(first_father), "/Exp_surr", name_phase, "/surr_");
                    
                    for (int s = 0; s < num_surr; s++)
                    {
                        addr[s + 1] = string.Concat(address_father, (s + 1).ToString());
                        Directory.CreateDirectory(addr[s + 1]);
                        
                    }
                    
                    for (int k = 0; k < number_electrodes; k++)
                    {   
                        //Console.WriteLine("numero de clusters: {0}",number_electrodes);
                        for (int s = 0; s < num_surr; s++)
                        {

                            string namef = namefile[k].Substring(namefile[k].LastIndexOf('/'));
                            namef = string.Concat(namef.Substring(0, namef.Length - 4), ".txt");
                            eee[s] = new FileStream(string.Concat(address_father, (s + 1).ToString(), namef), FileMode.Append);
                            sss[s] = new StreamWriter(eee[s]);

                            sss[s].WriteLine(n_samples);
                        }

                        FileStream electrode_k3 = new FileStream(string.Concat(namefile[k].Substring(0, namefile[k].Length - 4), ".bin"), FileMode.Create);
                        BinaryWriter El_k3 = new BinaryWriter(electrode_k3);
                        
                        
                        FileStream electrode_k2 = new FileStream(namefile[k], FileMode.Open);
                        StreamReader El_k2 = new StreamReader(electrode_k2);

                        //Console.WriteLine();

                        El_k2.ReadLine();
                        dim = int.Parse(El_k2.ReadLine());
                        int x;
                        
                        //El_k3.Write(n_samples);
                                                
                        for (int i = 0; i < num_surr; i++)
                        {
                            sss[i].WriteLine(dim);
                        }

                        //El_k3.Write(dim);
    
                        for (int i = 0; i < dim; i++)
                        {
                            x = int.Parse(El_k2.ReadLine());

                            //El_k3.Write(x);
                            
                            if (x - w > 0 && x + w < n_samples)

                                for (int surr = 0; surr < num_surr; surr++)
                                {
                                    sss[surr].WriteLine(rand.Next(x - w, x + w));
                                
                                }

                            else

                                for (int surr = 0; surr < num_surr; surr++)
                                {
                                    sss[surr].WriteLine(x);
                                }
                        }

                        for (int s = 0; s < num_surr; s++)
                        {
                            sss[s].Dispose();
                            eee[s].Close();
                        }

                        El_k2.Dispose();
                        electrode_k2.Close();
                        El_k3.Dispose();
                        electrode_k3.Close();
                        //File.Delete(namefile[k]);

                        

                    }
                    

                }


                #endregion
                

                #region if the files are binary

                else if (isbin == 1)

                {
                    Console.WriteLine("los datos son binarios");  
                }


                #endregion

                #endregion

                double count = 0;
                int bin_samples = (int)Math.Truncate(bin * fs / 1000);
                int number_bins = (int)Math.Truncate((double)n_samples / bin_samples);
                

                #region computation

                for (int id_k =0 ; id_k < id_electrodes.Length; id_k++)
                    filesDelete(address);

                double[,] mean_ncch = new double[2 * (int)L + 1,nchNS*(nchNS-1)/2];   
                double[,] devstd_ncch = new double[2 * (int)L + 1,nchNS*(nchNS-1)/2];   
                double[] mcc = new double[nchNS*(nchNS-1)/2];
                double[,] cc = new double[2 * (int)L + 1,nchNS*(nchNS-1)/2];
                double[ ] CC = new double[2 * (int)L + 1];
                double[,] mxc= new double[2 * (int)L + 1,nchNS*(nchNS-1)/2];
                double[,] devstd = new double[2 * (int)L + 1,nchNS*(nchNS-1)/2];
                double[,] TMC = new double[nchNS,nchNS];
                double vmedia =0;
                double[] Cross_matrix =new  double[number_electrodes*number_electrodes]; 

                double perc_tail_filter = 15.0d;//% de filtrado de cola de vetana definido como el 15%
                double perc_recompute = 40.0d; // % de ventana de recalculo

                int bin_for_filter = (int)( perc_tail_filter / 100.0 * L);
                int bin_for_re_compute = (int)( perc_recompute / 100.0 * L);

                int W = (int)L;

                string fecha_dat = Path.GetFileNameWithoutExtension((address).Substring(address.LastIndexOf("_") - 9));

                //FileStream xxoutmcc = new FileStream(string.Concat(path_mcc,"/","verificando_treshold_",fecha_dat,".txt"), FileMode.OpenOrCreate);
                //StreamWriter xxmcc_out = new StreamWriter(xxoutmcc);

                shufling_spdi objth = new shufling_spdi();
                shufling_spdi phip = new shufling_spdi();

                for (int i_col=0; i_col < nchNS*(nchNS-1)/2; i_col++)
                {                
                    for (int m_i=0 ; m_i < 2 * (int)L + 1; m_i++)
                    {

                        mean_ncch[m_i,i_col]=0;
                        devstd_ncch[m_i,i_col]=0; 
                           
                    }
                }

                Console.WriteLine("clusters activos: {0}",non_silent.Length);

                if(non_silent.Length!=0)
                {
                    for (int i=0 ; i < num_surr+1; i++)//System.Threading.Tasks.Parallel.For(0, num_surr + 1, i =>
                    {
                        
                        mxc = objth.correlationFunction(fsuavz,i,(int)L, mfr, window, silenti,bin_samples, addr[i] ,
                        nchNS, bin,fs,to_sort_electrodes,non_silent, n_thresh,path_delay,fecha_dat,
                        path_inhib, delay_inhib,recompute_inhib);

                        if(i==0)
                        {
                            cc =  mxc;         
                        }
                        else
                        { 
                            for (int i_col=0; i_col < nchNS*(nchNS-1)/2; i_col++)
                            {
                                for (int m_i=0 ; m_i < 2 * (int)L + 1; m_i++)
                                {   
                                    mean_ncch[m_i,i_col] +=  mxc[m_i,i_col];
                                    devstd_ncch[m_i,i_col] += Math.Pow(mxc[m_i,i_col],2); 
                                    
                                }
                                
                            }

                        }
                        
                    }

                    

                    double mean_CC;
                    
                    FileStream outmcc = new FileStream(string.Concat(path_mcc,"/","mcc_",fecha_dat,".txt"), FileMode.OpenOrCreate);
                    StreamWriter mcc_out = new StreamWriter(outmcc);

                    double[] mxtreshold_max =new  double[2 * (int)L + 1];
                    double[] mxtreshold_min =new  double[2 * (int)L + 1];
                    double [] zr =new double [2 * (int)L + 1];
                    if(num_surr > 2)
                    {
                        for (int i_col=0; i_col < nchNS*(nchNS-1)/2; i_col++)
                        {
                            for (int m_i=0 ; m_i < 2 * (int)L + 1; m_i++)
                            {

                                vmedia =  mean_ncch[m_i,i_col]/num_surr;
                                
                                double vespc = devstd_ncch[m_i,i_col]/num_surr;

                                devstd[m_i,i_col] = Math.Sqrt((vespc - Math.Pow(vmedia,2))*num_surr/(num_surr-1)); //para una muestra 
                            
                            }
                        }
                        
                        for (int i_col=0; i_col < nchNS*(nchNS-1)/2; i_col++)
                        {
                            for (int m_j=0; m_j < 2 * (int)L + 1; m_j++) 
                            {
                                vmedia =  mean_ncch[m_j,i_col]/num_surr;
        
                                double devstdr = devstd[m_j,i_col]; //desviacion standart 
            
                                CC[m_j]=  cc[m_j,i_col];

                                mxtreshold_max[m_j]=vmedia+2*devstdr;
                                mxtreshold_min[m_j]=vmedia-2*devstdr;
                            
                            }
                        
                            double max=0;
                            int index_max=0;

                            mean_CC = CC.Average();

                            double tempo = 0;

                            for (int o = 0; o < 2 * (int)L + 1; o++) //L = ampiezza finestra di correlazione in bin.
                            {
                                tempo = CC[o] - mean_CC;

                                if (Math.Abs(tempo) > Math.Abs(max))
                                {
                                    max = tempo;
                                    index_max = o;
                                }
                            }

                            if(recompute_inhib && max < 0 && Math.Abs(index_max - L) > W - bin_for_filter)
                            {
                                (max, index_max, mean_CC) = filter(max, index_max, (int)L, W, bin_for_filter, bin_for_re_compute, CC,  mean_CC);
                                
                            }

                        
                            //Treshold ok
                            if(max>0 && max > (mxtreshold_max[index_max]-mean_CC) && max > 0.001)                          
                                mcc[i_col]= max;//Math.Abs(max);

                            else if(max<0 && max < (mxtreshold_min[index_max]-mean_CC ))
                                mcc[i_col]=max;//Math.Abs(max);

                            else 
                                mcc[i_col]=0;
                                
                        }    
                    }
                    else
                    {

                        for (int i_col=0; i_col < nchNS*(nchNS-1)/2; i_col++)
                        {
                            for (int m_j=0; m_j < 2 * (int)L + 1; m_j++) 
                            {
                                CC[m_j]=  cc[m_j,i_col];  
                            }

                            double max=0;
                            int index_max=0;

                            mean_CC = CC.Average();
                            

                            double tempo = 0;

                            for (int o = 0; o < 2 * (int)L + 1; o++) //L = ampiezza finestra di correlazione in bin.
                            {
                                tempo = CC[o] - mean_CC;
                                
                                if (Math.Abs(tempo) > Math.Abs(max))
                                {
                                    max = tempo;
                                    index_max = o;
                                }
                            }

                            if(recompute_inhib && max < 0 && Math.Abs(index_max - L) > W - bin_for_filter)
                            {
                                (max, index_max,mean_CC) = filter(max, index_max, (int)L, W, bin_for_filter, bin_for_re_compute, CC,  mean_CC);
                                
                            }
                                        
                            mcc[i_col]= max;   
                            
                        } 

                    }
                    
                    int z,t,container=0;

                    for(int i = 0 ; i< nchNS; i++)
                    {
                        for(int p=0; p<i+1; p++)
                        {
                            if (i == p)
                                continue;

                            if (single_dimension[p] > single_dimension[i])
                            {
                                t = i;
                                z = p;
                            }
                            else
                            {
                                t = p;
                                z = i;
                            }
                            
                            Cross_matrix[non_silent[t] + non_silent[z] * number_electrodes] = mcc[container];
                            Cross_matrix[non_silent[z] + non_silent[t] * number_electrodes] = mcc[container];
                            container++; 
                        }
                    }
                    
                    for (int u = 0; u < number_electrodes - nchNS; u++)
                    {
                        for (int oo = 0; oo < number_electrodes; oo++)
                        {
                            Cross_matrix[silenti[u] * number_electrodes + oo] = 0;
                            Cross_matrix[oo * number_electrodes + silenti[u]] = 0;
                        }
                    }
                    
                    
                    for (int h = 0; h < number_electrodes; h++)
                    {
                        
                        for (int y = 0; y < number_electrodes; y++)
                        {
                            mcc_out.Write(string.Concat(Cross_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));
            
                        }
                        
                        mcc_out.WriteLine();

                    }

                    mcc_out.Close();

                   
                
                }
                //limpiar(address_parent,numero_de_experimentos,path_delay,path_inhib); 
                    #endregion

                    #region test
                    
                    #endregion

                
                // we completed an experiment 
                
              
            }

        
    }
    void test(ref float[] stats, double[] CM, string address_shuffl, int n)
    {


        string[] CM_shuffled = File.ReadAllLines(address_shuffl);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                //in input

                //generata dai surrogati
                if (CM[i + n * j] < Math.Abs(double.Parse(CM_shuffled[i].Split('\t')[j])))
                    stats[i + n * j]++;
            }
        }
    }

    public double[,] correlationFunction(bool fsuavz,int a_k,int L, float mfr, float window, int[] silenti, int bin_samples, 
    string address, int nchNS, float bin, float fs, LibraryCorrelations.add_on_function to_sort_electrodes, int[] non_silent, 
    double n_thresh,string path_delay,string fecha_dat, string path_inhib, bool delay_inhib,bool recompute_inhib)

    {
        double perc_tail_filter = 15.0d;//% de filtrado de cola de vetana definido como el 15%
        double perc_recompute = 40.0d; // % de ventana de recalculo

        int bin_for_filter = (int)( perc_tail_filter / 100.0 * L);
        int bin_for_re_compute = (int)( perc_recompute / 100.0 * L);

        double mean_CC = 0;
        double mean_CC2 = 0;


        int[] id_el={0};
        
        #region reading_of_surrogate

        System.Globalization.CultureInfo ci = new System.Globalization.CultureInfo("es-ES");
        System.Threading.Thread.CurrentThread.CurrentCulture = ci;

        string name = Path.GetFileName(address);

        string[] namefile2 = Directory.GetFiles(address);

        string[] namefile = new string[namefile2.Length];

        int number_electrodes = namefile.Length;

        // It's convenient to sort the name of the electrodes' file in crescent order, in order to semplify
        //the computing of the Connectivity Matrix (CM) and to better order the cross and partia correlogram

        int[] id_electrodes = new int[namefile.Length];

        // At first we save all the electrodes indexes in the order that they are red by the software

        for (int i = 0; i < number_electrodes; i++)
            id_electrodes[i] = int.Parse(Path.GetFileNameWithoutExtension((namefile2[i]).Substring(namefile2[i].LastIndexOf("_") + 1)));
            
        // Then we sort the array in crescent order, by using a built-in function

        Array.Sort(id_electrodes);

        // Now we use a function that simply returns the position of the electrodes file with a specific index (id)

        id_el = id_electrodes;

        for (int i = 0; i < namefile.Length; i++)

        //In this way we built a string vector with the electrodes sorted in crescent order, and we will use this vector in the rest
        //of the code
        {
            int file_position = to_sort_electrodes.electrode_find_file_position(namefile2, id_electrodes[i].ToString());

            if (file_position == -1)

                namefile[i] = namefile2[i];

            else

                namefile[i] = namefile2[file_position];

        }


        FileStream electrode = new FileStream(namefile[0], FileMode.Open);
        StreamReader El = new StreamReader(electrode);
        int n_files = namefile.Length;
        int n_samples = int.Parse(El.ReadLine());
        El.Dispose();
        electrode.Close();
        int temp = 0;
        int dim_total_spikes = 0;

        //Console.WriteLine(n_files);

        //array containing the index of last spike for every electrode in the cumulative array of spike times
        int[] dim_electrode2 = new int[number_electrodes];

        #endregion


        //We compute the correlation window's size

        double max;
        int index_max;

        #region reading files

        LibraryCorrelations.add_on_function funzione = new LibraryCorrelations.add_on_function();
        double normalization_factor;

        //array containing the number of spikes of each electrode
        int[] single_dimension = new int[number_electrodes];

        for (int k = 0; k < number_electrodes; k++)
        {
            FileStream electrode_k = new FileStream(namefile[k], FileMode.Open);
            StreamReader El_k = new StreamReader(electrode_k);
            // we determine the number of spikes (i.e. the number of elements contained in each file) 
            El_k.ReadLine();
            //the second element in every file is the number of spikes
            temp = int.Parse(El_k.ReadLine());

            //build non-silent and silent electrode indexes and build vector of single and cumulative dimensions

            if (temp / (n_samples / fs) >= mfr)
            {

                single_dimension[k] = temp;
                dim_electrode2[k] = dim_total_spikes + temp;
                //update the cumulative dimension
                dim_total_spikes = dim_total_spikes + temp;

            }


            El_k.Dispose();
            electrode_k.Close();
        }


        dim_electrode2 = funzione.nonzeros(dim_electrode2);
        single_dimension = funzione.nonzeros(single_dimension);

        int bin_spike = 0;

        //array containg the (occupied)bin times for the whole set of electrodes
        //i.e. if the bin of value -n indicates that a spike occured in the interval [n*(1/fs), (n+1)*(1/fs)]

        int[] spikes_data_2 = new int[dim_total_spikes];



        int[] dim_electrode = new int[nchNS];

        //array containg the index of the first spike for every electrode in the cumulative array of spikes
        int[] dim_electrode3 = new int[nchNS];

        for (int i = 1; i < nchNS; i++)
            dim_electrode3[i] = dim_electrode2[i - 1];

        #endregion

        int temp3 = 0;

        //filling the array of (occupied)bin times for the whole electreode's set
        for (int k = 0; k < nchNS; k++)
        {

            FileStream electrode_k2 = new FileStream(namefile[non_silent[k]], FileMode.Open);
            StreamReader El_k2 = new StreamReader(electrode_k2);

            El_k2.ReadLine();
            El_k2.ReadLine();
            bin_spike = (int.Parse(El_k2.ReadLine()) - 1) / (bin_samples);
            spikes_data_2[dim_electrode3[k]] = bin_spike;

            //(TODO: da controllare) ci manca assegnazione iniziale: 
            temp3 = bin_spike;
            //altrimenti: se uno spike successivo al primo ricadesse ancora nel primo bin, 
            //verrebbe inserito nel secondo bin perche con i = 1 vale sempre: bin_spike != temp

            for (int i = 1; i < single_dimension[k]; i++)
            {
                bin_spike = (int.Parse(El_k2.ReadLine()) - 1) / (bin_samples);

                //check if the spike time is inside the same time bin (we only consider if a spike occured in a bin)
                if (bin_spike != temp3)
                {
                    spikes_data_2[dim_electrode3[k] + i] = bin_spike;

                    temp3 = bin_spike;

                
                }

            }


            //number of spikes of last BIN for every electrode
            //number_spikes_per_bin[dim_electrode3[k] + single_dimension[k] - num] = num;

            int[] temp2 = funzione.nonzeros(spikes_data_2);

            //in single_dimension we store the size (in BIN) of each electrode k. Thus, we
            //compute it as temp2 (total number of spikes binned for electrodes from 0 to k) minus the 
            //same number for electrode from 0 to k-1 => number of binned spikes for electrode k.                        
            single_dimension[k] = temp2.Length - dim_electrode[k];

            //id dim_electrode we store the index of the first BIN of each electrode k.
            //Thus, the index of the first bin for k=0 is 0, the index of the first bin for k=1 is 0+single_dimension[0]
            if (k != nchNS - 1)
                dim_electrode[k + 1] = temp2.Length;

            El_k2.Dispose();
            electrode_k2.Close();

            
        }

        spikes_data_2 = funzione.nonzeros(spikes_data_2);

        dim_electrode[0] = 0;


        #region Output_preparation

        string crosscorrfolder = string.Concat(address, "/ncc_", mfr.ToString(), "-", bin.ToString(), "-", window.ToString(), "msec", "/", name).Replace(",", ".");

        Directory.CreateDirectory(crosscorrfolder);

        #endregion

        #region cross_time_computing

        //----------------------------------------------------------
        //----------------------------------------------------------
        // CROSS CORRELATION COMPUTING
        //----------------------------------------------------------              
        //----------------------------------------------------------


        string[] mxLhin = new string[number_electrodes*number_electrodes];
        
        double[,] mxCC = new double[2 * (int)L + 1,nchNS*(nchNS-1)/2];

        for (int mxl_inh=0; mxl_inh<number_electrodes*number_electrodes;mxl_inh++)
        {
            mxLhin[mxl_inh]=string.Concat(0);
        }

        double[] Cross_matrix = new double[number_electrodes * number_electrodes];
        int[] Delay_matrix = new int[number_electrodes * number_electrodes];
        int w = (int)(L);
        int[] idy = new int[nchNS];

        double[] CC = new double[2 * (int)L + 1];

        for (int i = 0; i < nchNS; i++)
        {
            idy[i] = int.Parse(Path.GetFileNameWithoutExtension((namefile[non_silent[i]]).Substring(namefile[non_silent[i]].LastIndexOf("_") + 1)));
        
        }

        int idy_max=idy.Max();

        int t, z;
        int spikes_x, spikes_y;
        int id_x, id_y;
        int flag;

        int r_k=0; 

        FileStream r_Table_id = new FileStream(string.Concat(path_inhib,"/id_",fecha_dat,"_inhib.txt"), FileMode.Append);
        StreamWriter r_Table = new StreamWriter(r_Table_id);

        for (int i = 0; i < nchNS; i++)
        {   
            //////////////////////////////////////////////////////////////////

            for (int p = 0; p < i + 1; p++)
            {
                flag = 0;

                if (i == p)
                    continue;

                if (single_dimension[p] > single_dimension[i])
                {
                    t = i;
                    z = p;
                }
                else
                {
                    t = p;
                    z = i;
                }

                id_x = dim_electrode[z];  //posizione di  partenza degli spikes di x
                id_y = dim_electrode[t];  //posizione di partenza degli spikes di y
                spikes_y = single_dimension[t]; //numero spikes di y
                spikes_x = single_dimension[z];
                //only to verify if the thread has been interrupted by the user
                

                //annulliamo Cxy
                for (int counter_for_Cxy = 0; counter_for_Cxy < 2 * (int)L + 1; counter_for_Cxy++)
                    CC[counter_for_Cxy] = 0;

                //scorriamo tutti i bin e ci posizioniamo in corrispondenza ad uno con almeno uno spike

            


                //scorrimento delgli spikes di y
                while (id_y < dim_electrode[t] + spikes_y || flag != 0)
                {
                    // esco se ho raggiunto l'ultimo spike di y e non ci sono più spike di y che cadono in finestre di correlazione di x (flag = 0);

                    if (flag != 0 && id_y == dim_electrode[t] + spikes_y)
                    {
                        // se invece ho finito gli spikes di y ma ce ne sono alcuni che cadono nella finestra di correlazione di spikes di x
                        id_y = flag;
                        flag = 0;
                        id_x++;
                        if (id_x == spikes_x + dim_electrode[z])
                        {
                            id_y = spikes_y + dim_electrode[t];
                            continue;
                        }
                    }
                    // controllo se lo spike idy di y cade dopo la finestra di correlazione dello spike idx di x
                    if (spikes_data_2[id_y] > spikes_data_2[id_x] + w)
                    {
                        //aggiustamento dell'indice idy dal quale riparte l'analisi in base al flag settato prima

                        if (flag != 0)
                        {
                            id_y = flag;
                            flag = 0;
                        }

                        id_x++; // analizzo il prossimo spike di reference x
                                // se gli spikes di reference sono finiti, esco mettendo idy = spikes_y +1 ed esco completamente;

                        if (id_x == spikes_x + dim_electrode[z])
                        {
                            id_y = spikes_y + dim_electrode[t];
                            continue;
                        }

                        continue;
                    }


                    //end if spike idy cade do po la finestra di correlazione di spike idx
                    //se lo spike idy non cade nella finestra dello spike idx passo ad analizzare lo spike successivo di y ed esco

                    if (spikes_data_2[id_y] < spikes_data_2[id_x] - w)
                    {
                        id_y++;
                        continue;
                    }


                    // Se lo spike id_y è nella finestra di correlazione, incremento il conteggio relativo al bin in questione espresso come y[id_y] – x[id_x]
                    CC[spikes_data_2[id_y] - spikes_data_2[id_x] + w]++;


                    // controllo che questo spike di y non cada nella finestra di correlazione dello spike successivo id_x + 1 di x
                    // se ciò si verifica, allora dovrò ripartire da questo spike di y al prossimo giro.

                    // se sto per uscire, e lo spike di x + w cade oltre il numero di bin, devo comunque analizzare un eventuale ultimo flag. 


                    if (flag != 0)
                    {
                        id_y++;
                        continue;
                    }

                    if (id_x + 1 == spikes_x + dim_electrode[z])
                    {
                        id_y++;
                        continue;
                    }

                    if ((spikes_data_2[id_y] >= (spikes_data_2[id_x + 1] - w)) && (spikes_data_2[id_y] <= (spikes_data_2[id_x + 1] + w)))
                    {
                        flag = id_y;
                    }

                    // analizzo lo spike successivo di y
                    id_y++;

                }  // end while

                normalization_factor = Math.Sqrt(spikes_x * spikes_y);
                max = 0;
                index_max = 0;

                for (int o = 0; o < 2 * (int)L + 1; o++) //L = ampiezza finestra di correlazione in bin.
                {
                    CC[o] = CC[o] / normalization_factor;
                }                
   

                /*
                se aplica un suavisado a la funcion de correlacion cruzada, a_k es igual a cero 
                se aplica el suavizado de lo contrario no se aplica.

                */

                //## inicia suavizado
                
                if(fsuavz && a_k==0)
                {
                    int nLeft = 10;
                    int nRight = 10;
                    int ngrp =3 ;

                    
                    DoubleVector NCC = new DoubleVector( CC.Length );

                    for ( int a_i = 0; a_i < CC.Length; a_i++ )
                    {
                        double a_CC = CC[a_i];
                        NCC[a_i] = a_CC;
                    }

                    DoubleVector vec = MovingWindowFilter.SavitzkyGolayCoefficients( nLeft, nRight, ngrp );
                    MovingWindowFilter filter = new MovingWindowFilter( nLeft, nRight, vec );
                    DoubleVector NCCH = filter.Filter(NCC, MovingWindowFilter.BoundaryOption.DoNotFilterBoundaryPoints);

                    for ( int a_ri = 0; a_ri < NCCH.Length; a_ri++ )
                    {
                        double a_NCC = NCCH[a_ri];
                        CC[a_ri] = a_NCC;
                    }

                }//## fin suavizado
                
                mean_CC = CC.Average();

                double tempo = 0;

                for (int o = 0; o < 2 * (int)L + 1; o++) //L = ampiezza finestra di correlazione in bin.
                {
                    tempo = CC[o] - mean_CC;
                    
                    if (Math.Abs(tempo) > Math.Abs(max))
                    {
                        max = tempo;
                        index_max = o;
                    }
                }
               
                if(recompute_inhib && max < 0 && Math.Abs(index_max - L) > w - bin_for_filter)
                {
                   (max, index_max,mean_CC) = filter(max, index_max, L, w, bin_for_filter, bin_for_re_compute, CC,  mean_CC);
                   
                }
                
                if(delay_inhib && max<0)
                {
                    r_Table.WriteLine("{0}{1}{2}",idy[p],".",idy[i]);
                }
                
                mxLhin[non_silent[t] + non_silent[z] * number_electrodes]=string.Concat(idy[p],".",idy[i]);

                Cross_matrix[non_silent[t] + non_silent[z] * number_electrodes] = max;//Math.Abs(max);
                //Console.WriteLine("{0}",non_silent[t] + non_silent[z] * number_electrodes);
                Delay_matrix[non_silent[t] + non_silent[z] * number_electrodes] = index_max - (int)L;
               
                mxLhin[non_silent[z] + non_silent[t] * number_electrodes]=string.Concat(idy[p],".",idy[i]);
                //Console.WriteLine("{0}\t{1}\t{2}",non_silent[z],  non_silent[t] , number_electrodes);
                Cross_matrix[non_silent[z] + non_silent[t] * number_electrodes] = max;//Math.Abs(max);
                Delay_matrix[non_silent[z] + non_silent[t] * number_electrodes] = -index_max + (int)L;

                //Console.WriteLine("{0}{1}{2}{3}{4}{5}",crosscorrfolder, "/r", name, '_', idy[p], ".txt");
                // falta clusteer ojo --> 486
                FileStream r_table = new FileStream(string.Concat(crosscorrfolder, "/r", name, '_', idy[p], ".txt"), FileMode.Append);
                FileStream rrr_table = new FileStream(string.Concat(crosscorrfolder, "/r", name, '_',idy_max, ".txt"), FileMode.Append);
                StreamWriter R_TABLE = new StreamWriter(r_table);

                //Saving the cross-correlation files
                R_TABLE.Write(idy[i]);
                R_TABLE.WriteLine();
                
                // Generate the coefficients.
                for (int q =0; q<2*(int)L+1;q++)
                {   
                    mxCC[q,r_k]=CC[q];   
                } 
               
                r_k++;

               

                for (int q = 0; q < 2 * (int)L + 1; q++)
                {
                   
                    R_TABLE.Write((CC[q]).ToString().Replace(",", "."));

                    R_TABLE.Write(R_TABLE.NewLine);
                }

                R_TABLE.Dispose();
                r_table.Close();
                
            }
        
        }
        
        r_Table.Dispose();
        r_Table_id.Close();


        for (int u = 0; u < number_electrodes - nchNS; u++)
        {
            for (int oo = 0; oo < number_electrodes; oo++)
            {

                mxLhin[silenti[u] * number_electrodes + oo]=string.Concat(0);
                mxLhin[oo * number_electrodes + silenti[u]]=string.Concat(0);

                Cross_matrix[silenti[u] * number_electrodes + oo] = 0;
                Cross_matrix[oo * number_electrodes + silenti[u]] = 0;
            }
        }

        FileStream output = new FileStream(string.Concat(crosscorrfolder, "_CC_Symmetric.txt"), FileMode.OpenOrCreate);
        FileStream output2 = new FileStream(string.Concat(crosscorrfolder, "_CC_Symmetric.bin"), FileMode.OpenOrCreate);
        FileStream output4 = new FileStream(string.Concat(crosscorrfolder, "_Delay.txt"), FileMode.OpenOrCreate);
        FileStream output3 = new FileStream(string.Concat(crosscorrfolder, "_CC_directed.txt"), FileMode.OpenOrCreate);
        FileStream output5 = new FileStream(string.Concat(crosscorrfolder, "_CC_directed.bin"), FileMode.OpenOrCreate);

    
        FileStream output6 = new FileStream(string.Concat(path_inhib,"/" ,"Links_inhib_",fecha_dat,".txt"), FileMode.OpenOrCreate);
        StreamWriter Out6 = new StreamWriter(output6);

        FileStream outdelay = new FileStream(string.Concat(path_delay,"/","Delay_",fecha_dat,".txt"), FileMode.OpenOrCreate);
        StreamWriter deyout = new StreamWriter(outdelay);
    


        StreamWriter Out = new StreamWriter(output);
        StreamWriter Out3 = new StreamWriter(output3);
        BinaryWriter Out2 = new BinaryWriter(output2);
        StreamWriter Out4 = new StreamWriter(output4); 
        BinaryWriter Out5 = new BinaryWriter(output5);
        


        //coincidence index
        //FileStream output3 = new FileStream(string.Concat(crosscorrfolder, "\\CC_CI.bin"), FileMode.OpenOrCreate);
        // BinaryWriter Out3 = new BinaryWriter(output3);
        
        //double[] MC_NCC = new double[number_electrodes*number_electrodes];      
        
        double[] CM_threshed = funzione.nonzeros(Cross_matrix);
        /*
        double threshold =MathNet.Numerics.Statistics.Statistics.Mean(CM_threshed) + n_thresh * MathNet.Numerics.Statistics.Statistics.StandardDeviation(CM_threshed);
        
        for (int g = 0; g < number_electrodes * number_electrodes; g++)
            if (Cross_matrix[g] < threshold)
                Cross_matrix[g] = 0;
        */        
        //int mxcont=0;

        for (int h = 0; h < number_electrodes; h++)
        {
            
            for (int y = 0; y < number_electrodes; y++)
            {
                Out.Write(string.Concat(Cross_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));

                //MC_NCC[mxcont] = Cross_matrix[h + y * number_electrodes]; 
                
                //mxcont++;

                Out2.Write(Cross_matrix[h + y * number_electrodes]);

                if (Delay_matrix[h + y * number_electrodes] > 0)

                    Cross_matrix[h + y * number_electrodes] = 0;

                Out4.Write(string.Concat(Delay_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));
                
                if(delay_inhib)
                {
                    Out6.Write(string.Concat(mxLhin[h + y * number_electrodes], '\t'));
                    deyout.Write(string.Concat(Delay_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));
                }
                Out3.Write(string.Concat(Cross_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));
                Out5.Write(Cross_matrix[h + y * number_electrodes]);
                
                
            }
            
            Out.WriteLine();
            Out3.WriteLine();
            Out4.WriteLine();

            if(delay_inhib)
            {
                deyout.WriteLine();
                Out6.WriteLine();
            }
        }

        Out.Close();
        Out2.Close();
        Out3.Close();
        Out4.Close();
        Out5.Close();

        if(delay_inhib)
        {
            deyout.Close();
            Out6.Close();
        }

        #endregion

        

        return mxCC;
    }


    public (double,int,double) filter(double max, int index_max,int L, int w,int bin_for_filter,int bin_for_re_compute,double [] CC, double mean_CC)
    { 
        double tempo=0; double  mean_CC2=0;
        
        bool real_peak = false;

        if (index_max - L < 0)
        {
            for (int k = 0; k <  bin_for_filter; k++)
                if (CC[k] - mean_CC > 0)
                {
                    real_peak = true;
                    break;
                }
        }
        else
        {
            for (int b = 2 * w  - bin_for_filter; b < 2 * w + 1; b++)
                if (CC[b] - mean_CC > 0)
                {
                    real_peak = true;
                    break;
                }
        }

        //recompute the peak in the eak recomputation window
        if (!real_peak)
        {   

            max = 0;
            mean_CC2 = 0;

            for (int o = L - bin_for_re_compute; o < L + bin_for_re_compute; o++) //L = ampiezza finestra di correlazione in bin.
            {
                mean_CC2 = CC[o] + mean_CC2;
            }
            mean_CC2 = mean_CC2 / (2 * (L - bin_for_re_compute));

            for (int o = L - bin_for_re_compute; o < L + bin_for_re_compute; o++) //L = ampiezza finestra di correlazione in bin.
            {
                tempo = CC[o] - mean_CC2;

                if (Math.Abs(tempo) > Math.Abs(max))
                {
                    max = tempo;
                    index_max = o;
                }
            }
        }
    
        return (max, index_max,mean_CC2);
    }




    public  double prueba_z(double x)
    {
        // constants
        double a1 =  0.254829592;
        double a2 = -0.284496736;
        double a3 =  1.421413741;
        double a4 = -1.453152027;
        double a5 =  1.061405429;
        double p  =  0.3275911;

        // Save the sign of x
        int sign = 1;
        if (x < 0)
            sign = -1;
        x = Math.Abs(x)/Math.Sqrt(2.0);

        // A&S formula 7.1.26
        double t = 1.0/(1.0 + p*x);
        double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.Exp(-x*x);

        return 0.5*(1.0 + sign*y);
    }


    }

}   
