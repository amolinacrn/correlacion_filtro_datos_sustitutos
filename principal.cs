
using System;
using System.Collections.Generic;
using System.Linq;
using System;
using System.IO;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using LibraryCorrelations;
//using  ClassicNCCH;
using ShuflingSpdi;
using FilteredFNCCH;


namespace  MainNCCH 
{
  	class main_NCCH 
  	{
    
		public main_NCCH (){}
         
       	static void Main() 
    	{   
            shufling_spdi delshuf = new shufling_spdi();  
            shufling_spdi dir_dt = new shufling_spdi();  

            string path_delay = @"MxDelays";
            string path_inhib = @"MxInhib";
            string [] set_exp = new string[]{@"Experiments"};//,@"ExpJan29",@"ExpDez20",@"ExpMar04",@"ExpJan14" ,@"ExpMar07"};

            for(int nNexp=0; nNexp < set_exp.Length; nNexp++)
            {
                string surr_path = set_exp[nNexp];
                
                string[] address_parent_principal = Directory.GetDirectories(surr_path);
                
                for(int i=0; i<address_parent_principal.Length;i++)
                { 
                    string[] address_parent_surr = Directory.GetDirectories(address_parent_principal[i]);

                    string address = address_parent_principal[i].ToString();
                    

                    string name_exp = address.Substring(address.LastIndexOf("_")+1);
                    
                    dir_dt.create_dir_dt(name_exp);
                    

                    int numero_de_experimentos =address_parent_surr.Length;
                    
                    bool wrdat=true;
                    bool delay_inhib=false;
                    bool inh_filtrar=false;
                    bool fsuavz=false;
                    float window=30f; //ms
                    int num_surr=50;
                    float fs=30000.0f;//(HZ)
                    float mfr=0.1f; //sṕikes/s
                    float bin=1.0f;//ms.
                    int wsurrgte=10;
                
                    //classic_NCCH cctime = new classic_NCCH();  
                    shufling_spdi shuf = new shufling_spdi();  
                    filtered_FNCCH fsurr  = new filtered_FNCCH();    
                    filtered_FNCCH fccn  = new filtered_FNCCH(); 
                
                    //cctime.computeCrossTime(address_parent, numero_de_experimentos,fs,mfr,bin,window);

                    shuf.shufling(fsuavz,address_parent_surr, numero_de_experimentos, fs,  mfr,  bin,  window,  num_surr, wsurrgte,name_exp,path_delay,path_inhib,delay_inhib, inh_filtrar);
                    //fsurr.filtar_surr_correlations(address_parent_surr,numero_de_experimentos, num_surr, bin, wsurrgte,wrdat);
                    //fccn.filtar_cross_correlation(address_parent_surr,numero_de_experimentos, num_surr, bin, wsurrgte,wrdat);
                   
                }
                
                 for(int i=0; i<address_parent_principal.Length;i++)
                 { 
                     string[] address_parent_surr = Directory.GetDirectories(address_parent_principal[i]);

                     int numero_de_experimentos =address_parent_surr.Length;

                     delshuf.limpiar(address_parent_surr,numero_de_experimentos,path_delay,path_inhib);
                }
                 
            }
        }
        
    }
}


