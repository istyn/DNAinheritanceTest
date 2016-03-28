/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Name                    :Isaac Styles
// Department Name : Computer and Information Sciences 
// File Name              :Program.cs
// Purpose                 :determine a strand of best match for each chromosome,
//						        determine a parent/child relation, and phase out the missing parent.
//							
// Author			: Isaac Styles, styles@goldmail.etsu.edu
// Create Date	: Dec. 3, 2014
//
//--------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Modified Date	: Dec 1, 2015
// Modified By		: Isaac Styles
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Project3
{
	class Program
	{
		static void Main(string[] args)
		{
			//int[,]  _strands ={{1,1},{1,2},{2,1},{2,2}};
			int[] d1ChromoPointers, d2ChromoPointers;                                         //pointer locations to first allele of each chromosome
			//alorithm counts chromosome length during file parsing and ensures there are the same number of compares
			//USED TO guarantee that chromosomes are of same length prior to comparison,
			string[,] d1 = ReadDNAFile(@"J:\DNAB.txt", out d1ChromoPointers);     //string[,] columns are chromosomeNumber, strand 1, strand 2
			string[,] d2 = ReadDNAFile(@"J:\DNAM.txt", out d2ChromoPointers);
			string diagTextOut = "";                                        //contains main diagnostic output before initializing output file
			double[,] strandScores = new double[25, 4];                     //holds the scores of all 4 comparisons of each chromosome
			double[,] chromoScore = new double[25,2];                       //holds the resulting product of each persons probability of inheritance 
			double[] parentScore = {1d, 1d};                                //holds the percent match for person 1 and person 2
			double correlation;                                             //the ratio of correlation of relation
			char[] missingParent;                                           //holds the array of chromoNumber and Allele of the extrapolated parent
			int[] pMissingParentChromo = new int[25];                       //1-based pointer to non-matched strand of child
			int pPointerToParent = -1;                                      //0-based pointer to person determined to be parent
			bool d1Female = false;
			bool d2Female = false;

			for (int i = 0; i < 26; i++)                                    //check that chromosomes are equilength
			{
				if (d1ChromoPointers[i] != d2ChromoPointers[i])
				{
					throw new Exception("Chromosomes not aligned at chromo " + (i + 1));
				}
			}
			for (int c = 0; c < 25; c++)        //compare one strand to the other persons strands, looking for possibility of inheritance
			{
				for (int comparisonIndex = 0; comparisonIndex < 4; comparisonIndex++) //do 4 comparisions for each chromosome
				{
					int match = 0, total = 0;

					for (int i = 0; i < d1ChromoPointers[c + 1] - d1ChromoPointers[c]; i++)
					{
						if (comparisonIndex == 0)                       //person 1 strand 1
						{
							if (d1[i + d1ChromoPointers[c], 1] == d2[i + d1ChromoPointers[c], 1] || d1[i + d1ChromoPointers[c], 1] == d2[i + d1ChromoPointers[c], 2])
							{
								match++;
							}
							total++;
						}
						else if (comparisonIndex == 1)                  //person 1 strand 2
						{
							if (d1[i + d1ChromoPointers[c], 2] == d2[i + d1ChromoPointers[c], 1] || d1[i + d1ChromoPointers[c], 2] == d2[i + d1ChromoPointers[c], 2])
							{
								match++;
							}
							total++;
						}
						else if (comparisonIndex == 2)                  //person 2 strand 1
						{
							if (d2[i + d1ChromoPointers[c], 1] == d1[i + d1ChromoPointers[c], 1] || d2[i + d1ChromoPointers[c], 1] == d1[i + d1ChromoPointers[c], 2])
							{
								match++;
							}
							total++;
						}
						else if (comparisonIndex == 3)                  //person 2 strand 2
						{
							if (d2[i + d1ChromoPointers[c], 2] == d1[i + d1ChromoPointers[c], 1] || d2[i + d1ChromoPointers[c], 2] == d1[i + d1ChromoPointers[c], 2])
							{
								match++;
							}
							total++;
						}
						else { throw new Exception("Invalid comparison pointer. Valid values are: 0, 1, 2, 3"); }

					}                    
					strandScores[c, comparisonIndex] = match / (double)total;
				}
			}
			for (int i = 0; i < 25; i++)                                //multiply probability of inheritance for each persons' pair of chromosomes
			{
				chromoScore[i, 0] = strandScores[i, 0] * strandScores[i, 1];
				chromoScore[i, 1] = strandScores[i, 2] * strandScores[i, 3];
			}
			for (int i = 0; i < 25; i++)                                //multiply probability of inheritance for each persons' genome
			{
				parentScore[0] = parentScore[0] * chromoScore[i, 0];
				parentScore[1] = parentScore[1] * chromoScore[i, 1];
			}
			//calculation percent different of inheritance scores
			correlation = Math.Abs(parentScore[0]-parentScore[1]);
			correlation = correlation /((parentScore[0] + parentScore[1]) / 2);
			string parent = "";                                             //temp string for parent name
			string child = "";                                              //temp for child name
			if (parentScore[0] < parentScore[1])                            //display the results
			{
				parent = "Person 1";
				child = "Person 2";
				pPointerToParent = 0;

			}
			else if (parentScore[0] > parentScore[1])
			{
				parent = "Person 2";
				child = "Person 1";
				pPointerToParent = 1;
			}
			else                                                            //scored the same
			{
				throw new Exception("Best-fit for parent/child relation not found. Check that two unique individuals are being compared.");
			}
					/*Display results in debug window*/
			diagTextOut += "      Parent / Child Relationship:\r\n  Parent:" + parent + "\r\n   Child:" + child + "\r\n";
			diagTextOut += "         Correlation (larger is better): " + correlation;
			diagTextOut += "\r\n              s1+t1  s1+t2  s2+t1  s2+t2";    //display chromosome scores
			for (int i = 0; i < 25; i++)
			{
				diagTextOut += "\r\nChromosome" + (i+1) + ":";
				for (int j = 0; j < 4; j++)
				{
					diagTextOut += Math.Round(strandScores[i, j], 4).ToString().PadLeft(7, ' ');

				}
				diagTextOut += "\r\n";
			}
			System.Diagnostics.Debug.WriteLine(diagTextOut);
			diagTextOut = "";

					/*Find pointers to best matched chromosome of child*/
			if (pPointerToParent == 1)                          //person 2 is parent
			{
				for (int i = 0; i < 25; i++)
				{
					if (strandScores[i,0]>strandScores[i,1])      //if P1 chromo 1 best match
					{
						pMissingParentChromo[i] = 2;                    //chromo 2 is from missing parent
					}
					else                                        //if P1 chromo 2 best match
					{
						pMissingParentChromo[i] = 1;                    //chromo 1 is from missing parent
					}                    
				} 
			}
			else                                                //person 1 is parent
			{
				for (int i = 0; i < 25; i++)                
				{
					if (strandScores[i, 2] > strandScores[i, 3])  //P2 chromo 1 best match
					{
						pMissingParentChromo[i] = 2;                    //chromo 2 is best match
					}
					else                                        //P2 chromo 2 is best match
					{
						pMissingParentChromo[i] = 1;                    //chromo 1 is best match
					}
				} 
			}
																/*Begin output of missing parent*/
			missingParent = new char[d1ChromoPointers[d1ChromoPointers.GetUpperBound(0)]]; //initialize the missing parent 
																/*copy child's non-matched strand to missing parent*/
			if (pPointerToParent == 1)                          //person 2 is parent
			{
				for (int i = 0; i < 23; i++)
				{
					for (int j = 0; j < d1ChromoPointers[i + 1] - d1ChromoPointers[i]; j++)
					{
						missingParent[j + d1ChromoPointers[i]] = d1[j + d1ChromoPointers[i], pMissingParentChromo[i]][0];
					}
				} 
			}
			else
			{
				for (int i = 0; i < 23; i++)
				{
					for (int j = 0; j < d1ChromoPointers[i + 1] - d1ChromoPointers[i]; j++)
					{
						missingParent[j + d1ChromoPointers[i]] = d2[j + d1ChromoPointers[i], pMissingParentChromo[i]][0];
					}
				} 
			}
			/* Determine sex of parent and child, then output sex chromosomes*/
			d1Female = IsFemale(ref d1ChromoPointers, ref d1);
			d2Female = IsFemale(ref d2ChromoPointers, ref d2);

			if (pPointerToParent == 1)                          //person 2 is parent
			{
				if (d2Female)     //parent is female; give missing parent X and look for Y
				{
					if (d1Female) //both parent and child are female, NO KNOWN Y chromosome
					{
						diagTextOut = "Exporting Missing Parent: No known Y chromosome exists. Defaulting to '0'.";
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)//copy child's Y to missing parent
						{
							missingParent[i] = '0';
						}
					}
					else                                        //person 2 is female parent, but Y exists in P1
					{
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)//copy child's Y to missing parent
						{
							missingParent[i] = d1[i, pMissingParentChromo[i]][0];
						}
					}
				}
				else                                            //person 2 is male parent; give missing parent X & no Y
				{
					if (d1Female) //child is female; give child's chromosome 24
					{
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)
						{
							missingParent[i] = d1[i, pMissingParentChromo[i]][0];
						}
					}
					else                                        //both parent and child are male, NO KNOWN !Y chromosome
					{
						diagTextOut = "Exporting Missing Parent: No known !Y chromosome exists. Defaulting to '0'.";
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)
						{
							missingParent[i] = '0';
						}
					}
				}
			}
			else                                                //Person 1 is parent
			{
				if (d1Female)     //parent is female; look for Y
				{
					if (d2Female) //both parent and child are female, NO KNOWN Y chromosome
					{
						diagTextOut = "Exporting Missing Parent: No known Y chromosome exists. Defaulting to '0'.";
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)//copy child's Y to missing parent
						{
							missingParent[i] = '0';
						}
					}
					else                                        //person 1 is female parent, but Y exists in P2
					{
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)//copy child's Y to missing parent
						{
							missingParent[i] = d2[i, pMissingParentChromo[23]][0];
						}
					}
				}
				else                                            //person 2 is male parent; give missing parent X & no Y
				{
					if (d2Female) //child is female; give child's chromosome 24
					{
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)
						{
							missingParent[i] = d2[i, pMissingParentChromo[23]][0];
						}
					}
					else                                        //both parent and child are male, NO KNOWN !Y chromosome
					{
						diagTextOut = "Exporting Missing Parent: No known !Y chromosome exists. Defaulting to '0'.";
						for (int i = d1ChromoPointers[23]; i < d1ChromoPointers[24]; i++)
						{
							missingParent[i] = '0';
						}
					}
				}
			}
																/*copy mitochondrial dna*/
			if (pPointerToParent == 1)                          //person 2 is parent
			{
				for (int i = d1ChromoPointers[24]; i < d1ChromoPointers[25]; i++)
				{
					missingParent[i] = d1[i, pMissingParentChromo[24]][0];
				} 
			}
			else                                                //person 1 is parent
			{
				for (int i = d1ChromoPointers[24]; i < d1ChromoPointers[25]; i++)
				{
					missingParent[i] = d2[i, pMissingParentChromo[25]][0];
				} 
			}
			System.Diagnostics.Debug.WriteLine(diagTextOut);    //write the diagnostic output to the debug window
			for (int i = 0; i <= missingParent.GetUpperBound(0); i++)//output the best guess for missing parent
			{
				Console.WriteLine(missingParent[i].ToString());
			}
		}

		/// <summary>
		/// Determines whether the given DNA is female.
		/// </summary>
		/// <param name="chromoPointers">The chromo pointers.</param>
		/// <param name="dna">The dna.</param>
		/// <returns>TRUE if female. FALSE if male.</returns>
		private static bool IsFemale(ref int[] chromoPointers, ref string[,] dna)
		{
			int count = 0;
			int total = 0;

			for (int i = chromoPointers[23]; i < chromoPointers[24]; i++)
			{
				if ("0" == dna[i, 1])                        //chromosome 24 is male Y
				{
					count++;
				}               
				total++;
			}
			if (count / (double)total > .9d)                //if majority of Y was 0
			{
				return true;
			}
			else return false;
		}

		/// <summary>
		/// Determines the Y haplogroup.
		/// </summary>
		/// <param name="dna">The dna.</param>
		/// <param name="chromoPointers">The chromo pointers.</param>
		/// <returns></returns>
		public static List<string> DetermineHaplogroup(ref string[,] dna, ref int[] chromoPointers)
		{
			List<string> haplogroup = new List<string>();       //list of haplogroups that were expressed in Y chromosome
			string[] h;                                         //list of most specific haplogroups
			bool[] ignore;                                      //corresponds to found haplogroups that may be ignored
			string[,] haplogroups = ReadHaplotypeFile(@"../../haplogroups.txt");    //haplogroups to search for
			for (int i = chromoPointers[23]; i < chromoPointers[24]; i++)             //look at Y chromosome
			{
				for (int j = 0; j < haplogroups.GetUpperBound(0); j++)  //iterate through the haplogroups
				{
					if (dna[i, 3] == haplogroups[j, 2])                  //position matches known haplogroup
					{
						if (dna[i, 1] == haplogroups[j, 3])              //mutation matches haplogroup
						{
							if (!haplogroup.Contains(haplogroups[j, 0]))
							{
								haplogroup.Add(haplogroups[j, 0]);
							}
						}
					}
				}
			}
			h = haplogroup.ToArray();                           //put haplogroups into array to maintain an index
			ignore = new bool[h.Length];
			for (int x = 0; x < h.Length; x++)                  //compare each pair of haplogroups once
			{
				for (int y = x + 1; y < h.Length; y++)
				{
					if (h[x] == null || h[y] == null || h[x].Length == h[y].Length)  //if haplogroup is less specific, or haplogroups have same number of characters
					{
						continue;                               //then cannot determine specificity
					}
					if (h[x].Length < h[y].Length)              //if lengths differ, begin looking at characters to see if they belong to same branch
					{
						bool morespecific = true;               //set to false when a character doesn't match, indicating two different branches
						for (int i = 0; i < h[x].Length; i++)   //iterate the chars, looking for a mismatch
						{
							if (h[x][i] != h[y][i])
							{
								morespecific = false;           //if mismatch is found, then the longer haplogroup is not more specific
							}

						}
						if (morespecific == true)               //if mismatch is not found, then the longer haplogroup is more specific
						{
							h[x] = null;                   //ignore the shorter haplogroup
						}
					}
					else                                        //the second haplogroup is longer, compare chars and look for a mismatch
					{
						bool morespecific = true;
						for (int i = 0; i < h[y].Length; i++)
						{
							if (h[x][i] != h[y][i])
							{
								morespecific = false;
							}

						}
						if (morespecific == true)
						{
							h[y] = null;
						}
					}
				}
			}
			return haplogroup;
		}
		public static string[,] ReadHaplotypeFile(string filePath)
		{
			string[,] output = new string[1333, 4];    //preallocated array for Haplogroups
			string theLine = "";
			int x = 0;
			try
			{
				using (StreamReader sr = new StreamReader(filePath))
				{
					string[] fields;
					while (theLine != null)     //infinite loop with counter x
					{
						theLine = sr.ReadLine();
						if (theLine != null && theLine[0] != '#' && theLine[2] != 'i')
						{
							fields = theLine.Split(new char[] { '\t' });
							output[x, 0] = fields[0];               //group
							output[x, 1] = fields[1];               //rsid
							output[x, 2] = fields[2];               //position
							output[x, 3] = fields[3];               //mutation

							x++;
						}
					}
				}
			}
			catch (IOException e)
			{
				Console.WriteLine("The file could not be read: ");
				Console.WriteLine(e.Message);
			}

			return output;
		}

		/// <summary>
		/// Reads the dna file, returning a table with columns={ ChromosomeNumber, Strand 1, and Strand 2 }
		/// Comments begin line with #.
		/// </summary>
		/// <param name="filePath">The file path of DNA.</param>
		/// <param name="chromoPointers">Returns array of pointers to each chromosome.</param>
		/// <returns>table with columns={ ChromosomeNumber, Strand 1, and Strand 2 }</returns>
		public static string[,] ReadDNAFile(string filePath, out int[] chromoPointers)
		{
			chromoPointers = new int[26];                      //pointers to beginning of each chromosome (25+ endPointer)
			int x = 0, i=0;                             //pointer in string[], and int[] chromoP
			string[,] output = new string[701478,4];    //preallocated array for ChrisDNA.txt
			string theLine = "";

			try
			{
				using (StreamReader sr = new StreamReader(filePath))
				{
					string[] fields;
					while (theLine != null)     //infinite loop with counter x
					{
						theLine = sr.ReadLine();
						if (theLine!=null && theLine[0] != '#' && theLine[2] != 'i')
						{
							fields = theLine.Split(new char[] { '\t' });
							output[x,0] = fields[1];                //chromosome Number
							output[x, 3] = fields[2];               //rsid AT THE END
							output[x, 1] = fields[3];               //Strand 1
							output[x, 2] = fields[4];               //strand 2
													/*This if statment saves a pointer to the beginning of every chromosome.
														Allows comparison of chromosome lengths. */
							if (x == 0) 
							{
								chromoPointers[0] = 0;
							}
							else if (output[x, 0] != output[x - 1, 0])
								{
									chromoPointers[++i] = x;
								}
							x++;
						}
					}
					chromoPointers[chromoPointers.GetUpperBound(0)]=x;

				}


			}
			catch (IOException e)
			{
				Console.WriteLine("The file could not be read: ");
				Console.WriteLine(e.Message);
			}
			
			return output;
		}

	}
}
