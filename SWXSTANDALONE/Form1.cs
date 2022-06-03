using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Diagnostics;
using System.IO;
using System.Threading;

namespace SWXSTANDALONE
{
    public partial class Form1 : Form
    {
        public int file_index = 0;
        public Form1()
        {
            //this.TopMost = true;
            InitializeComponent();
            
        }
        
        private void button1_Click(object sender, EventArgs e)
        {

            file_index++;
            label2.Text = "Extracting Weldments!";
            label2.ForeColor = System.Drawing.Color.Red;
            SolidWorksSingleton.Holla(file_index);
            label2.Text = "Weldments extracted";
            label2.ForeColor = System.Drawing.Color.Green;
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            String selection = comboBox1.GetItemText(this.comboBox1.SelectedItem);
            if (selection.Contains("UR"))
            {
                pictureBox1.Load(@"C:\Users\Matei\source\repos\SWXSTANDALONE\SWXSTANDALONE\ur1.png");
            }
            else
            {
                pictureBox1.Load(@"C:\Users\Matei\source\repos\SWXSTANDALONE\SWXSTANDALONE\moto1.png");
            }
            

        }

        private void button2_Click(object sender, EventArgs e)
        {
            //ProcessStartInfo start = new ProcessStartInfo();
            //start.FileName = @"C:\Users\Matei\\anaconda3\python.exe";
            //start.Arguments = @"C:\Users\Matei\OneDrive\Desktop\weaving\weaving3.py";
            //start.UseShellExecute = false;
            //start.RedirectStandardOutput = true;
            //using (Process process = Process.Start(start))
            //{
            //    using (StreamReader reader = process.StandardOutput)
            //    {
            //        string result = reader.ReadToEnd();
            //        Debug.Print(result);
            //    }
            //}
            

            new Thread(() =>
            {
                ProcessStartInfo start = new ProcessStartInfo();
                start.FileName = @"C:\RoboDK\bin\RoboDK.exe";
                start.Arguments = @"C:\Users\simat\OneDrive\Desktop\ElementaryOperations\YaskawaPositionerThingDemo.rdk";
                start.UseShellExecute = false;
                start.RedirectStandardOutput = true;
                using (Process process = Process.Start(start))
                {
                    using (StreamReader reader = process.StandardOutput)
                    {
                        string result = reader.ReadToEnd();
                        Debug.Print(result);
                    }
                }
            }).Start();

            Thread.Sleep(2000);
            new Thread(() =>
            {
                ProcessStartInfo start = new ProcessStartInfo();
                start.FileName = @"C:\Users\simat\.conda\envs\EWO\pythonw.exe";
                start.Arguments = @"C:\Users\simat\OneDrive\Desktop\ElementaryOperations\Horn\hrn5.py";
                start.UseShellExecute = false;
                start.RedirectStandardOutput = true;
                using (Process process = Process.Start(start))
                {
                    using (StreamReader reader = process.StandardOutput)
                    {
                        string result = reader.ReadToEnd();
                        Debug.Print(result);
                    }
                }
            }).Start();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            new Thread(() =>
            {
                ProcessStartInfo start = new ProcessStartInfo();
                start.FileName = @"C:\Users\simat\.conda\envs\EWO\pythonw.exe";
                start.Arguments = @"C:\Users\simat\OneDrive\Desktop\ElementaryOperations\Horn\tst3.py";
                start.UseShellExecute = false;
                start.RedirectStandardOutput = true;
                using (Process process = Process.Start(start))
                {
                    using (StreamReader reader = process.StandardOutput)
                    {
                        string result = reader.ReadToEnd();
                        Debug.Print(result);
                    }
                }
            }).Start();



            //String selection = comboBox1.GetItemText(this.comboBox1.SelectedItem);
            //if (selection.Contains("UR"))
            //{
            //    new Thread(() =>
            //    {
            //        ProcessStartInfo start = new ProcessStartInfo();
            //        start.FileName = @"C:\RoboDK\bin\RoboDK.exe";
            //        start.Arguments = @"C:\Users\Matei\OneDrive\Desktop\ElementaryOperations\URStationBucket.rdk";
            //        start.UseShellExecute = false;
            //        start.RedirectStandardOutput = true;
            //        using (Process process = Process.Start(start))
            //        {
            //            using (StreamReader reader = process.StandardOutput)
            //            {
            //                string result = reader.ReadToEnd();
            //                Debug.Print(result);
            //            }
            //        }
            //    }).Start();

            //    Thread.Sleep(3000);
            //    void message() {
            //        string message = "Unable to generate trajectories. Consider changing the robot.";
            //        string caption = "Error";
            //        MessageBoxButtons buttons = MessageBoxButtons.OK;
            //        DialogResult result;

            //        result = MessageBox.Show(message, caption, buttons);
            //        if (result == System.Windows.Forms.DialogResult.OK)
            //        {
            //            this.Close();
            //        }
            //    }
            //    message();

            //}
            //else
            //{
            //    new Thread(() =>
            //    {
            //        ProcessStartInfo start = new ProcessStartInfo();
            //        start.FileName = @"C:\Users\Matei\\anaconda3\python.exe";
            //        start.Arguments = @"C:\Users\Matei\OneDrive\Desktop\ElementaryOperations\executesimulation.py";
            //        start.UseShellExecute = false;
            //        start.RedirectStandardOutput = true;
            //        using (Process process = Process.Start(start))
            //        {
            //            using (StreamReader reader = process.StandardOutput)
            //            {
            //                string result = reader.ReadToEnd();
            //                Debug.Print(result);
            //            }
            //        }
            //    }).Start();
            //}


        }

        private void button4_Click(object sender, EventArgs e)
        {

            new Thread(() =>
            {
                ProcessStartInfo start = new ProcessStartInfo();
                start.FileName = @"C:\Users\simat\.conda\envs\EWO\pythonw.exe";
                start.Arguments = @"C:\Users\simat\OneDrive\Desktop\ElementaryOperations\Horn\tst4.py";
                start.UseShellExecute = false;
                start.RedirectStandardOutput = true;
                using (Process process = Process.Start(start))
                {
                    using (StreamReader reader = process.StandardOutput)
                    {
                        string result = reader.ReadToEnd();
                        Debug.Print(result);
                    }
                }
            }).Start();
            //new Thread(() =>
            //{
            //    ProcessStartInfo start = new ProcessStartInfo();
            //    start.FileName = @"C:\Users\Matei\\anaconda3\python.exe";
            //    start.Arguments = @"C:\Users\Matei\OneDrive\Desktop\ElementaryOperations\executeprogram.py";
            //    start.UseShellExecute = false;
            //    start.RedirectStandardOutput = true;
            //    using (Process process = Process.Start(start))
            //    {
            //        using (StreamReader reader = process.StandardOutput)
            //        {
            //            string result = reader.ReadToEnd();
            //            Debug.Print(result);
            //        }
            //    }
            //}).Start();

        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {

        }
    }
}
