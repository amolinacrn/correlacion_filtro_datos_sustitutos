using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.IO;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using LibraryCorrelations;

namespace  ShuflingApplicationIsi
{
  	class ShuflingSurrogIsi
  	{
    
		public ShuflingSurrogIsi(){}
        /* 
       	static void Main()
    	{
           for(int i=0; i<2;i++)
           {
            string[] DirName = new string[]{@"Exp_sua_2020Mar4/Exp_surr",@"Exp_sua_2020Mar4/Exp_sua_2020Mar4/NCC_0.1-2-25msec"};
            // Checking if Directory Exist or not
            if (Directory.Exists(DirName[i]))
            {
 
                // This will delete the
                // Directory if it is empty
                Directory.Delete(DirName[i],true);
 
            }
           }
            DirectoryInfo di = new DirectoryInfo(@"Exp_sua_2020Mar4");
            FileInfo[] files = di.GetFiles();

            foreach (FileInfo file in files)
            {
                file.Delete();
            }
           
			ShuflingSurrog computshufling = new ShuflingSurrog();
            computshufling.shufling();  
            
		}*/

        void filesDelete(int id_clust)
        {

            string addresdir = @"Exp_sua_2020Mar4/Exp_sua_2020Mar4";

            string FileName = string.Concat(addresdir,"/sua_2020Mar4_",id_clust,".bin");
            //Console.WriteLine(FileName);
            // Checking file exists or not
            if (File.Exists(FileName))
            {
 
            // Deleting the file
            File.Delete(FileName);
  
            // Checking if the file is deleted
            // successfully or not
                  
            }
        } 

      
        private void shufling()
        {
            
            //#region Initialization

            string address = "o";
            int[] id_el={0};
			int numero_de_experimentos =1;
			int isbin=0;
		
			float fs=30000.0f; // (Hz) 
			float mfr=0.1f; //(spikes/s) 
			float bin=2f; //(ms) 
			float window=25.0f;//(ms)	
            double n_thresh =0.0d;
            int num_surr = 5;                   
            int wsurrgte = 25;
            //float alfa = 0.05f;


			//donde, fs: creo que es la frecuencia de muestreo, window: es el ancho de la ventana, bin: es el ancho de bin
			// mfr = spikes/s=0.1:frecuencia de disparo.
		
			string[] address_parent = new string[] {@"Exp_sua_2020Mar4/Exp_sua_2020Mar4"}; //grafmejorajs"};
			int number_of_experimental_set = address_parent.Length;
		
                

                for (int num_set = 0; num_set < numero_de_experimentos; num_set++)
                {

                    address = address_parent[num_set].ToString();

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
                   
                        int w = (int)(wsurrgte * (fs / 1000));
                        //addr[0] = address;
                        var rand = new Random();
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



                            El_k2.ReadLine();
                            dim = int.Parse(El_k2.ReadLine());
                            int[] x = new int[dim] ;
                            int[] ISI = new int[dim-1] ;

                            //El_k3.Write(n_samples);

                                                       

                            for (int i = 0; i < num_surr; i++)
                            {
                                sss[i].WriteLine(dim-1);
                            }

                            //El_k3.Write(dim);

                            for (int i = 0; i < dim; i++)
                                x[i] = int.Parse(El_k2.ReadLine());                                    
                                           
                            for (int i_n = 0; i_n < dim-1; i_n++)
                                ISI[i_n] = x[i_n+1]-x[i_n];                                                      
                            
                          
                            for (int surr = 0; surr < num_surr; surr++)
                            {
                                var randomized = ISI.OrderBy(item => rand.Next());

                                foreach(var nxsurr in randomized)
                                {   
                                    int mIndex = rand.Next(ISI.Length);
                                    
                                    sss[surr].WriteLine(nxsurr);
                                    
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
                        filesDelete(id_electrodes[id_k]);    
                                    
                                    
                    System.Threading.Tasks.Parallel.For(0, num_surr + 1, i =>
                    {
                        correlationFunction( (int)L, mfr, window, silenti,bin_samples, addr[i] ,nchNS, bin,fs,to_sort_electrodes,non_silent, n_thresh);
                    });
                   
                 
                    #endregion

                    #region test
                  
                    #endregion

                   
                    // we completed an experiment 
                    
                  
                }
            

            string[] address_log = new string[] {address};

            #region log_file_writing

            string log_folder = address_log[0] ;
            Directory.CreateDirectory(log_folder);
            FileStream LOG = new FileStream(string.Concat(log_folder, "Log_PC_CC.txt"), FileMode.Create);
            StreamWriter log = new StreamWriter(LOG);

            log.WriteLine("Parameters used in the correlation analysis");
            log.WriteLine();
            log.WriteLine(string.Concat("Bin size (ms)", '\t', '\t', '\t', '\t', bin.ToString()));
            log.WriteLine(string.Concat("Correlation Window (ms)", '\t', '\t', '\t', window.ToString()));
            log.WriteLine(string.Concat("Sampling Frequency (Hz)", '\t', '\t', '\t', fs.ToString()));
            log.WriteLine(string.Concat("minimum mean firing rate (spike/s)", '\t', mfr.ToString()));
            log.WriteLine(string.Concat("peak_train_folder", '\t', '\t', '\t', address));
            log.WriteLine(string.Concat("Electrodes Selected", '\n'));
            for (int log_i = 0; log_i < id_el.Length; log_i++)
            {
                log.WriteLine(string.Concat(id_el[log_i].ToString(), '\n'));
            }
            log.Dispose();
            LOG.Close();

            #endregion
     
            
                          
        }
            
        #region test
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

        #endregion

        void correlationFunction(int L, float mfr, float window, int[] silenti, int bin_samples, string address, int nchNS, float bin, float fs, LibraryCorrelations.add_on_function to_sort_electrodes, int[] non_silent, double n_thresh)
        {
            int[] id_el={0};
           
            #region reading_of_surrogate

            System.Globalization.CultureInfo ci = new System.Globalization.CultureInfo("es-ES");
            System.Threading.Thread.CurrentThread.CurrentCulture = ci;

            string name = Path.GetFileName(address);

            string[] namefile2 = Directory.GetFiles(address);

            string[] namefile = new string[namefile2.Length];

            int number_electrodes = namefile.Length;

            int[] n_id_clusters=new int[number_electrodes];

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

            string crosscorrfolder = string.Concat(address, "/NCC_", mfr.ToString(), "-", bin.ToString(), "-", window.ToString(), "msec", "/", name).Replace(",", ".");

            Directory.CreateDirectory(crosscorrfolder);

            #endregion

            #region cross_time_computing

            //----------------------------------------------------------
            //----------------------------------------------------------
            // CROSS CORRELATION COMPUTING
            //----------------------------------------------------------              
            //----------------------------------------------------------



            double[] Cross_matrix = new double[number_electrodes * number_electrodes];
            int[] Delay_matrix = new int[number_electrodes * number_electrodes];
            int w = (int)(L);
            int[] idy = new int[nchNS];

            double[] CC = new double[2 * (int)L + 1];

            for (int i = 0; i < nchNS; i++)
            {
                idy[i] = int.Parse(Path.GetFileNameWithoutExtension((namefile[non_silent[i]]).Substring(namefile[non_silent[i]].LastIndexOf("_") + 1)));
            
            }

            int t, z;
            int spikes_x, spikes_y;
            int id_x, id_y;
            int flag;

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

                    for (int o = 0; o < 2 * (int)L + 1; o++) //L = ampiezza finestra di correlazione in bin.
                    {               
                        if (CC[o]> max)
                        {
                            max = CC[o];
                            index_max = o;
                        }
                    }


                    Cross_matrix[non_silent[t] + non_silent[z] * number_electrodes] = max;
                    Delay_matrix[non_silent[t] + non_silent[z] * number_electrodes] = index_max - (int)L;

                    Cross_matrix[non_silent[z] + non_silent[t] * number_electrodes] = max;
                    Delay_matrix[non_silent[z] + non_silent[t] * number_electrodes] = -index_max + (int)L;


                  
                    FileStream r_table = new FileStream(string.Concat(crosscorrfolder, "/r", name, '_', idy[p], ".txt"), FileMode.Append);
                    StreamWriter R_TABLE = new StreamWriter(r_table);

                    //Saving the cross-correlation files
                    R_TABLE.Write(idy[i]);
                    R_TABLE.WriteLine();

                    for (int q = 0; q < 2 * (int)L + 1; q++)
                    {
                        R_TABLE.Write((CC[q]).ToString().Replace(",", "."));
                        //Console.WriteLine((CC[q]).ToString().Replace(",", "."));
                        //to remove


                        R_TABLE.Write(R_TABLE.NewLine);
                    }

                    R_TABLE.Dispose();
                    r_table.Close();
                    
                }
            
                      n_id_clusters[i]=idy[i];
            }

            for (int u = 0; u < number_electrodes - nchNS; u++)
            {
                for (int oo = 0; oo < number_electrodes; oo++)
                {

                    Cross_matrix[silenti[u] * number_electrodes + oo] = 0;
                    Cross_matrix[oo * number_electrodes + silenti[u]] = 0;
                }
            }

            FileStream output = new FileStream(string.Concat(crosscorrfolder, "_CC_Symmetric.txt"), FileMode.OpenOrCreate);
            FileStream output2 = new FileStream(string.Concat(crosscorrfolder, "_CC_Symmetric.bin"), FileMode.OpenOrCreate);
            FileStream output4 = new FileStream(string.Concat(crosscorrfolder, "_Delay.txt"), FileMode.OpenOrCreate);
            FileStream output3 = new FileStream(string.Concat(crosscorrfolder, "_CC_directed.txt"), FileMode.OpenOrCreate);
            FileStream output5 = new FileStream(string.Concat(crosscorrfolder, "_CC_directed.bin"), FileMode.OpenOrCreate);

            StreamWriter Out = new StreamWriter(output);
            StreamWriter Out3 = new StreamWriter(output3);
            BinaryWriter Out2 = new BinaryWriter(output2);
            StreamWriter Out4 = new StreamWriter(output4);
            BinaryWriter Out5 = new BinaryWriter(output5);

            //coincidence index
            //FileStream output3 = new FileStream(string.Concat(crosscorrfolder, "\\CC_CI.bin"), FileMode.OpenOrCreate);
            // BinaryWriter Out3 = new BinaryWriter(output3);

          
            double[] CM_threshed = funzione.nonzeros(Cross_matrix);
            double threshold =0;// MathNet.Numerics.Statistics.Statistics.Mean(CM_threshed) + n_thresh * MathNet.Numerics.Statistics.Statistics.StandardDeviation(CM_threshed);
            //Console.WriteLine(threshold);
            for (int g = 0; g < number_electrodes * number_electrodes; g++)
                if (Cross_matrix[g] < threshold)
                    Cross_matrix[g] = 0;
                    

            for (int h = 0; h < number_electrodes; h++)
            {
                for (int y = 0; y < number_electrodes; y++)
                {
                    Out.Write(string.Concat(Cross_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));
                    Out2.Write(Cross_matrix[h + y * number_electrodes]);

                    if (Delay_matrix[h + y * number_electrodes] > 0)

                        Cross_matrix[h + y * number_electrodes] = 0;

                    Out4.Write(string.Concat(Delay_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));
                    Out3.Write(string.Concat(Cross_matrix[h + y * number_electrodes].ToString(), '\t').Replace(",", "."));
                    Out5.Write(Cross_matrix[h + y * number_electrodes]);
                }

                Out.WriteLine();
                Out3.WriteLine();
                 Out4.WriteLine();
            }

            Out.Close();
            Out2.Close();
            Out3.Close();
            Out4.Close();
            Out5.Close();


            #endregion


        }


    }

}   