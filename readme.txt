
                         MatLab ECG Processing Demo


I. How to use:

   1. Copy files 

         ecgdemo.m;
         ecgdemowinmax.m;
         ecgdemodata1.mat;
         ecgdemodata2.mat

      (all the files except this "readme.txt" file) into MatLab's work
      directory.

   2. Start up MatLab.

   3. Type in 

         >> ecgdemo

      and press "Enter".

   4. You will get two figures everyone of which has 6 plots:

         a. The first plot shows original ECG data;
         b. The second plot contains corrected ECG - the low-frequency
            component was removed;
         c. The third plot shows the data after 1-st filtering pass, the
            filter window is of default size so the result is not "clear";
         d. The fourth plot shows detected peaks - on this stage some peaks
            can be skipped;
         e. Here we analyse the result of peak detection and optimize the
            filter window size. So the fifth plot contains the result of
            2-d filtering pass;
         f. Sixth plot shows the final result.

II. Package content:

   1. ecgdemo.m        - MatLab demo script;
   2. ecgdemowinmax.m  - MatLab window filter script;
   3. ecgdemodata1.mat - 1-st data sample;
   4. ecgdemodata2.mat - 2-d data sample;
   5. readme.txt       - this file.



                                Downloaded from:
                                   www.librow.com

                                Feedback:
                                   Sergey Chernenko: S.Chernenko@librow.com

    Surya modified the code by adding Heart Rate calculation  and some 
    additional comments. Contact: p.surya1994@gmail.com  

    Thanks p.surya1994@gmail.com for giving a meaningful project about ECG developement
    The file is updated with more ideal from our groups (Luong, Luan, Phong)
    We have reworked with files for more effiency in team work and changes for better understanding in general
    It includes about detect medical condition of persion related to R peak 
      + First, It is Atrial Fibrillation
      + We also update more data for testing
      + There are more added aglorithm to detect other cases of R peak.

