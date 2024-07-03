using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.IO;
using LibraryCorrelations;
using System.Drawing.Imaging;
using System;


namespace FilteredFNCCH
{
    class filtered_FNCCH
    {

        public filtered_FNCCH (){}


        public void filtar_cross_correlation(string[] xdir_surr,int numero_de_experimentos,int ndsurr, float xbin,int wid,bool wrdat)
    	{   
            for (int num_set = 0; num_set < numero_de_experimentos; num_set++)
            {
                //directorio donde se van a escribir todos los clusters
                string addir = xdir_surr[num_set].ToString();
               
                string fecha_dat = Path.GetFileNameWithoutExtension((addir).Substring(addir.LastIndexOf("_") - 9));
                
                string path = string.Concat(addir,"/cc",fecha_dat,"w",wid);//
               
                try
                {
                    // Determine whether the directory exists.
                    if (Directory.Exists(path))
                        Directory.Delete(path,true);
                    
                    // Try to create the directory.
                    if(wrdat)
                    {
                        DirectoryInfo di = Directory.CreateDirectory(path);
                        Console.WriteLine("{0}",path);
                    }
                    string[] xmfile = Directory.GetFiles(addir,"id_*"); 
                    
                    DirectoryInfo  deletfile = new DirectoryInfo(addir);    
                    
                    FileInfo[] files = deletfile.GetFiles();
                    
                    if(xmfile.Length !=0)
                    {    
                        foreach (FileInfo file in files)
                        {
                            file.Delete();
                        }
                    }
                    
                }

                catch (Exception e)
                {
                    Console.WriteLine("El proceso falló: {0}", e.ToString());
                }
                finally {}

               

                filtered_FNCCH importdatancch =new filtered_FNCCH();

                
        
                string crosfolder = string.Concat(path);
                
                string[] experimento = Directory.GetDirectories(addir,"Exp_sua*");
                      
                string[] dir_ccn = Directory.GetDirectories(experimento[0],"ncc*");   

                string[] xddir = Directory.GetDirectories(dir_ccn[0],"Exp_sua*");                    
                
                string FileName = string.Concat(xddir[0]);
                
               importdatancch.DataImport_NCCH(xbin, FileName,path, addir,crosfolder, wrdat);  
                
            
                
            }
		}
        
  	    public void filtar_surr_correlations(string[] xdir_surr,int numero_de_experimentos,int ndsurr, float xbin,int wid,bool wrdat)
    	{   
            for (int num_set = 0; num_set < numero_de_experimentos; num_set++)
            {
                //Console.WriteLine(numero_de_experimentos);
                //directorio donde se van a escribir todos los clusters
                string addir = xdir_surr[num_set].ToString();
                //Console.WriteLine(addir );
                string fecha_dat = Path.GetFileNameWithoutExtension((addir).Substring(addir.LastIndexOf("_") -9));
                
                string path = string.Concat(addir,"/sr",fecha_dat,"w",wid);//
                

                try
                {
                    // Determine whether the directory exists.
                    if (Directory.Exists(path))
                        Directory.Delete(path,true);
                    
                    // Try to create the directory.
                    if(wrdat)
                    {
                        DirectoryInfo di = Directory.CreateDirectory(path);
                        Console.WriteLine("{0}",path);
                    }
                    
                    string[] xmfile = Directory.GetFiles(addir,"id_*"); 
                    
                    DirectoryInfo  deletfile = new DirectoryInfo(addir);    
                    
                    FileInfo[] files = deletfile.GetFiles();
                    
                    if(xmfile.Length !=0)
                    {    
                        foreach (FileInfo file in files)
                        {
                            file.Delete();
                        }
                    }
                    
                }

                catch (Exception e)
                {
                    Console.WriteLine("El proceso falló: {0}", e.ToString());
                }
                finally {}


                filtered_FNCCH importdatancch =new filtered_FNCCH();

                
                for (int i=0; i < ndsurr; i++)
                {   
                    string crosfolder = string.Concat(path,"/surr_",i+1);
                                 
                    //directorio de donde queremos leer los datos
                    string[] addres = Directory.GetDirectories(addir,"*surr");
                    
                    string[] experimento = Directory.GetDirectories(addres[0],"Exp*");
                                    
                    //esta linea es para datos surrogate
                    string crosscorrfolder = string.Concat(experimento[0],"/surr_",i+1);
                    
                    string[] dir_ccn = Directory.GetDirectories(crosscorrfolder,"ncc*");                    
                    
                    string FileName = string.Concat(dir_ccn[0],"/surr_",i+1);
                    
                    importdatancch.DataImport_NCCH(xbin, FileName,path, addir,crosfolder, wrdat);  
                   
                
                }
            }
		}
    
        private void FNCCH_Function (string addres_ncch_id,string idcexp,double [] t, double [] x, int el1, int el2)
        {          
            bool filtrar=true;  
                           
            double W = (x.Length - 1)/2.0;
            int w = (int) Math.Round(W);

            FileStream r_Table_id = new FileStream(string.Concat(addres_ncch_id,"/id_",idcexp,"_inhib.txt"), FileMode.Append);
            StreamWriter r_Table = new StreamWriter(r_Table_id);
              
            //filtering tail effect eanbled for the "filtered cross-correlogram"
            //checked the filtered cross-correlogram analisys

            if( filtrar == true)
            {   

                
                //filtered cross-correlogram's peak
                double mean_x = x.Average(); //Console.WriteLine(mean_x);
                double filtered_max = x.Max() - mean_x;
                double filtered_min = x.Min() - mean_x;

                // aqui se calcula la media de la matriz NCCH (que es x) para calcular el enlace ihibitorio, asi: x.Average();
                //si el valor absoluto de (x.Min() - mean_x) es mayor al valor abosoluto de (x.Max() - mean_x), entonces, 
                //aqui hay un enlace inhibitorio y no un exitatorio como fue detectado por NCCH.
                // 


                //check if it we have a peak related to a max or a min value              
                //for inhibitory connections
                
                if (Math.Abs(filtered_min) > Math.Abs(filtered_max))
                {
                    //Console.WriteLine("holis cara de verga");

                    bool need_to_recompute = true;

                    double perc_tail_filter = 15.0d;//% de filtrado de cola de vetana definido como el 15%
                    double perc_recompute = 40.0d; // % de ventana de recalculo

                    //adding the mean subtracted (filtered) correlogram

                    int tail_x = (int)(2 * w - w * perc_tail_filter / 100);
                    int tail2_x = (int)(w * perc_tail_filter / 100);
                    int tail_y = (int)(w - w * perc_recompute / 100);
                    int tail2_y = (int)(w + w * perc_recompute / 100);
                    int bin_for_filter = (int) (w * perc_tail_filter / 100);

                    //searching the index of the negative peak identified in the "full" correlogram
                    int index_min = 0;
                    double val1 = 0;

                    for (int i = 0; i < x.Length; i++)
                    {
                        val1 = x[i] - mean_x;
                        //retreiving the index of the negative peak
                        if (val1 == filtered_min)
                            index_min = i;
                    }

                    //checking if the negative peak falls in the tail effect area
                    //if it is true the check if it is a FP or a "real" TP negative peak due to inhibitory connection
                    //the peak falls outside of the window

                    //de aqui para abajo es para verificar si el un pico cae en la de cola, si un pico cae en esta region
                    // se lo evalua si correponde a un enlace inhibitorio o no,
                    if (index_min >= tail2_x && index_min < tail_x)
                    {
                        need_to_recompute = false;
                        //Console.WriteLine("{0}\t{1}",el1,el2);
                        r_Table.WriteLine("{0}\t{1}",el1,el2);
                    }
                    else if (index_min >= 0 && index_min < tail2_x)
                    {
                        /////
                        for (int i = 0; i < bin_for_filter; i++) //(in tail effect window)
                        {
                            //if we find at least one positive value in the tail effect window it means that it is a real negative peak
                            if (x[i] - mean_x > 0)
                            {
                                need_to_recompute = false;
                                //Console.WriteLine("{0}\t{1}",el1,el2);
                                r_Table.WriteLine("{0}\t{1}",el1,el2);
                                break;
                               
                            }
                        }
                    }
                    else if(index_min >= tail_x && index_min < (2 * w + 1))
                    {
                        /////
                        for (int i = (2 * w + 1) - bin_for_filter; i < (2 * w + 1); i++)
                        {
                            //if we find at least one positive value in the tail effect window it means that it is a real negative peak
                            if (x[i] - mean_x > 0)
                            {
                                need_to_recompute = false;
                                //Console.WriteLine("{0}\t{1}",el1,el2);
                                r_Table.WriteLine("{0}\t{1}",el1,el2);
                                break;
                            }
                        }
                    }
                 
                    string len = "";
                    if ((x[index_min] - mean_x).ToString().Length > 7)
                        len = (x[index_min] - mean_x).ToString().Substring(0, 7);
                    else
                        len = (x[index_min] - mean_x).ToString();

                    //esta parte del codigo es para mirar si existe o no un enlace inhibitorio o no,
                    //se debe recalucular en una ventana mas pequeña para mirar si en la ventana de corrlacion
                    //mas pequeña se encuentra un enlace, que puede ser exitatorio o inhibitorio.
                    
                    //we need to recompute the peak in the peak recompute window                                     
                    if (need_to_recompute)
                    {
                        //computing the peak on the peak recompute area
                        double peak = 0;
                        double val = 0;
                        int idx_peak = 0;
                        double new_mean = 0;

                        //computing the new mean of the peak recompute window
                        for (int i = tail_y; i < tail2_y + 1; i++)
                        {
                            new_mean = x[i] + new_mean;
                        }
                        new_mean = new_mean / (tail2_y + 1 - tail_y);

                        //plot of the peak recompute window as a distinct series
                        for (int i = tail_y; i < tail2_y + 1; i++)
                        {
                            val = x[i] - new_mean;
                            
                            //effective peak's index recompute needet for the plot
                            if (Math.Abs(val) >= Math.Abs(peak))
                            {
                                peak = val;
                                idx_peak = i;
                            }

                        }
                        //new color related to the peak recompute window
                       

                        string conn = " ";
                        if(x[idx_peak] - new_mean > 0)
                        {
                            conn = "Excitatorio";
                            //Console.WriteLine("{0}:\t{1} --> {2}",conn,el1,el2);
                        }
                        else
                        {
                            conn = "Inhibitorio";
                            //Console.WriteLine("{0}\t{1}",el1,el2);
                            r_Table.WriteLine("{0}\t{1}",el1,el2);
                        }
                        if ((x[idx_peak] - new_mean).ToString().Length > 7)
                            len = (x[idx_peak] - new_mean).ToString().Substring(0, 7);
                        else
                            len = (x[idx_peak] - new_mean).ToString();

                        
                    }

                    // NO need to recompute the negative peak
                    else
                    {
                        //changing the text of the legend for the inhibitory peak
                        //string text = CHART.Series[number_Series + 4].LegendText;
                        //CHART.Series[number_Series + 4].LegendText = string.Concat(text, "\nNo Peek Recomputation");                
                    }

                }//end the peak is a Min

                else
                {
                    //CHART.Series.Add(string.Concat(el1.ToString(), " vs ", el2.ToString(), " NCCH (Cross-Correlogram)"));
                    //number_Series++;
                    //CHART.Series[number_Series].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
                    //for (int i = 0; i < x.Length; i++)
                    //    CHART.Series[number_Series].Points.AddXY(t[i], x[i]);

                    //CHART.Series.Add(string.Concat(el1.ToString(), " vs ", el2.ToString(), "  FCCH"));
                    
                    //the index of the positive peak of the filtered cross-correlogram
                    int index_max = 0;

                    double val = 0;
                    double max_val = 0;
                    for (int i = 0; i < x.Length; i++)
                    {
                        val = x[i] - mean_x;
                        //retreiving the index of the positive peak
                        if (val == filtered_max)
                        {
                            index_max = i;
                            max_val = val;
                        }
                            
                        
                    }
                    

                    //adding the positive valued peak of the filtered cross-correlogram related to Excitatory connection
                    string len = "";
                    if ((x[index_max] - mean_x).ToString().Length > 7)
                        len = (x[index_max] - mean_x).ToString().Substring(0, 7);
                    else
                        len = (x[index_max] - mean_x).ToString();
                }//end case of excitatory connection

                //r_Table.Write("{0}\t{1}",0,0);

               
            }//end filtered cross-correlogram   
                    
            else
            {
                //standard cross-correlogram's peak
                double standard_max = 0; //no identification of the inhibitory connections
                int index_max = 0;
               
                double val = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    val = x[i];
                    if  (val > standard_max)
                    {
                        standard_max = val;
                        index_max = i;
                    }                       
                }

                string len = "";
                if (x[index_max].ToString().Length > 7)
                    len = x[index_max].ToString().Substring(0, 7);
                else
                    len = x[index_max].ToString();
                //adding the positive valued peak related to Excitatory connection
            }   
            
            r_Table.Dispose();
            r_Table_id.Close();     
             
        }

        public void DataImport_NCCH(float ancho_bin, string addres,string addres_ncch, string addres_ncch_id,string crosscorrfolder,bool wrdat)
        {

            int el1 = 0;
            int el2 = 0;

            string[] namefile = Directory.GetFiles(addres);
            
            int[] id_clusters =new int[namefile.Length];
             
            for(int b_i = 0; b_i < namefile.Length; b_i++)
            {
                id_clusters[b_i] = int.Parse(Path.GetFileNameWithoutExtension((namefile[b_i]).Substring(namefile[b_i].LastIndexOf("_") + 1)));
               
            }
            
            int ord = 0;

            for (int a_k = 0; a_k < namefile.Length;a_k++)
            {
                for(int b_k = a_k; b_k < namefile.Length; b_k++)
                {
                 
                    if(id_clusters[b_k] < id_clusters[a_k])
                    {   
                        ord = id_clusters[a_k];
                        id_clusters[a_k] = id_clusters[b_k];
                        id_clusters[b_k] = ord;
                    }
             
                }
            }
            
            string idcexp = Path.GetFileNameWithoutExtension((addres_ncch_id).Substring(addres_ncch_id.LastIndexOf("_") - 9));
            string xdir_xadres=string.Concat(addres_ncch_id,"/id_",idcexp,".txt");
            
            FileStream r_table_id = new FileStream(xdir_xadres, FileMode.Append);
            StreamWriter R_TABLLE = new StreamWriter(r_table_id);

           for (int a_i = 0; a_i < namefile.Length ; a_i++ )
           {          
                for (int a_j = a_i+1; a_j < namefile.Length; a_j++)    
                {                      
                    
                    try
                    { 
                        el1 = id_clusters[a_i];
                        el2 = id_clusters[a_j];
                    }
                    catch
                    {
                        Console.WriteLine("Debe introducir el índice entero de los electrodos a trazar", "Error de digitación de los electrodos");
                        return;
                    }

                    if(el1 == el2)
                    {
                        Console.WriteLine("No hay gráfico de autocorrelograma", "Error de digitación de los electrodos");
                        return;
                    }

                    int min_el = Math.Min(el1, el2);
                    int max_el = Math.Max(el1, el2);

                    int n = namefile.Length;
                    int[] id_electrodes = new int[namefile.Length];

                    for (int i = 0; i < n; i++)
                    {
                        try
                        {
                        id_electrodes[i] = int.Parse(Path.GetFileNameWithoutExtension((namefile[i]).Substring(namefile[i].LastIndexOf("_") + 1)));
                        }               
                        catch
                        {
                        //if there are other fles in the same folder (e.g., 'CC_Symmetric.bin')
                        id_electrodes[i] = 0;
                        continue;                   
                        }
                    }
                
                    add_on_function var10 = new add_on_function();

                    string pre_File = "";

                    try
                    {
                        pre_File = namefile[var10.find(id_electrodes, min_el)];
                    }
                    catch
                    {
                        Console.WriteLine(string.Concat("Electrodo ", min_el.ToString(), "no encontrado"), "Error de selección de electrodos");
                        return;
                    }
                    
                    //first electrode's index
                    FileStream electrodefile = new FileStream(pre_File, FileMode.Open);
                    StreamReader contenuto = new StreamReader(electrodefile);
                    
                    //starting position for the correlogram related to electrode el1 -> el2
                    int pos = 1;
                    //total samples dimension of the correlogram
                    int dim = 0;
                    double val = 0;
                   val = double.Parse(contenuto.ReadLine());

                    while (val != max_el)
                    {
                        pos++;
                        try
                        {
                            val = double.Parse(contenuto.ReadLine());
                        
                        }
                        catch
                        {
                            Console.WriteLine(string.Concat("El correlograma no tiene muestras (", max_el.ToString() ,"es probablemente una silenciosa)"), "Error de selección de electrodos");
                            contenuto.Dispose();
                            electrodefile.Close();
                            return;
                        }
                
                    }

                    //counting the number of samples of the correlogram
                    val = double.Parse(contenuto.ReadLine());
                    
                    //while reachiang a integer value or the value we read is 0
                    while (val == 0 || val != Math.Truncate(val))
                    {
                        dim++;
                        try
                        {
                            val = double.Parse(contenuto.ReadLine());
                            
                        }
                        catch
                        {
                            //end of file, it can be the case of the correlogram for the last electrode
                            break;
                        }   
                         
                    }   
                         
                    //initializing the correlogram vector  
                    double[] x = new double [dim];

                    electrodefile.Close();
                    contenuto.Dispose();

                    FileStream electrodefile2 = new FileStream(pre_File, FileMode.Open);
                    StreamReader contenuto2 = new StreamReader(electrodefile2);
                    
                    //positioning to the electrodes's el2 start index
                    for (int z = 0; z < pos; z++)
                    {
                        contenuto2.ReadLine();
                        

                    }
                    //reading el2 data
                    for (int i = 0; i < dim; i++)
                    {
                        x[i] = double.Parse(contenuto2.ReadLine());
                       
                    }

                    electrodefile2.Close();
                    contenuto2.Dispose();

                    double [] t = new double[dim];
                    float bin_size;           
                    try
                    {
                        bin_size = ancho_bin;
                    }

                    catch
                    {
                        Console.WriteLine("Insertar el tamaño del bin", "Bin_Size_Insert");
                        return;
                    }

                    //filling the X-axis: time (ms)

                    for (int k = 0; k < dim; k++)
                    {
                        t[k] = (-(dim - 1) / 2.0 + k)*bin_size;
                  
                    } 

                    FNCCH_Function(addres_ncch_id,idcexp,t, x, el1, el2);
                          
                    if(wrdat)
                    {
                        Directory.CreateDirectory(crosscorrfolder);

                        FileStream r_table = new FileStream(string.Concat(crosscorrfolder,"/rsurr_",el1,"_",el2,".txt"), FileMode.Append);
                        StreamWriter R_TABLE = new StreamWriter(r_table);
                    
                        for (int q = 0; q < dim; q++)
                        {
                            R_TABLE.Write("{0}\t\t{1}",t[q].ToString(),(x[q]).ToString());
                             
                            R_TABLE.Write(R_TABLE.NewLine);
                        }

                        R_TABLE.Dispose();
                        r_table.Close();

                        R_TABLLE.WriteLine("{0}\t{1}",el1,el2);
                    }
                }
            }

            R_TABLLE.Write("{0}\t{1}",0,0);

            R_TABLLE.Dispose();
            r_table_id.Close();
       
        }
        
    }
}
