// --------------------------------------------------------------------------------------------------------------------
// <copyright file="VonMisesMaterial3D.cs" company="National Technical University of Athens">
//   To be decided
// </copyright>
// <summary>
//   Class for 3D Von Mises materials.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

using System;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Continuum
{
	/// <summary>
	///   Class for 3D Von Mises materials.
	/// </summary>
	/// <a href = "http://en.wikipedia.org/wiki/Von_Mises_yield_criterion">Wikipedia -Von Mises yield criterion</a>
	public class VonMisesMaterial3D : IIsotropicContinuumMaterial3D
	{
		private const string PLASTIC_STRAIN = "Plastic strain";
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";
		private GenericConstitutiveLawState currentState;

		/// <summary>
		///   The Poisson ratio value of an incompressible solid.
		/// </summary>
		private const double PoissonRatioForIncompressibleSolid = 0.5;

		/// <summary>
		///   The total number of strains.
		/// </summary>
		private const int TotalStrains = 6;

		/// <summary>
		///   The total number of stresses.
		/// </summary>
		private const int TotalStresses = TotalStrains;

		/// <summary>
		///   An array needed for the formulation of the consistent constitutive matrix.
		/// </summary>
		private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrix = new[,]
			{
				{  2.0/3.0, -1.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, 2.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, -1.0/3.0, 2.0/3.0, 0,   0,   0  },
				{  0,  0,  0, 1.0, 0,   0   },
				{  0,  0,  0, 0,   1.0, 0   },
				{  0,  0,  0, 0,   0,   1.0 }
			};

		/// <summary>
		///   The constitutive matrix of the material while still in the elastic region.
		/// </summary>
		private readonly Matrix elasticConstitutiveMatrix;

		/// <summary>
		///   Hardening modulus for linear hardening.
		/// </summary>
		private readonly double hardeningModulus;

		/// <summary>
		///   The Poisson ratio.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
		/// </remarks>
		private double poissonRatio;

		/// <summary>
		///   The shear modulus.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Shear_modulus">Wikipedia - Shear Modulus</a>
		/// </remarks>
		private readonly double shearModulus;

		/// <summary>
		///   The yields stress.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Yield_%28engineering%29">Yield (engineering)</a>
		///   The yield strength or yield point of a material is defined in engineering and materials science as the stress at which a material begins to deform plastically.
		/// </remarks>
		private readonly double yieldStress;

		/// <summary>
		///   The young modulus.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
		/// </remarks>
		private double youngModulus;

		/// <summary>
		///   The constitutive matrix of the material.
		/// </summary>
		private Matrix constitutiveMatrix;

		/// <summary>
		///   The array of incremental strains.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
		/// </remarks>
		private double[] incrementalStrains = new double[6];

		/// <summary>
		///   Indicates whether this <see cref = "IStructuralMaterial" /> is modified.
		/// </summary>
		private bool modified;

		/// <summary>
		///   The current plastic strain.
		/// </summary>
		private double plasticStrain;

		/// <summary>
		///   The new plastic strain.
		/// </summary>
		private double plasticStrainNew;

		/// <summary>
		///   The array of stresses.
		/// </summary>
		private double[] stresses = new double[6];

		/// <summary>
		///   The array of new stresses.
		/// </summary>
		private double[] stressesNew = new double[6];

		/// <summary>
		///   Initializes a new instance of the <see cref = "VonMisesMaterial3D" /> class.
		/// </summary>
		/// <param name = "youngModulus">
		///   The young modulus.
		/// </param>
		/// <param name = "poissonRatio">
		///   The Poisson ratio.
		/// </param>
		/// <param name = "yieldStress">
		///   The yield stress.
		/// </param>
		/// <param name = "hardeningModulus">
		///   The hardening ratio.
		/// </param>
		/// <exception cref = "ArgumentException"> When Poisson ratio is equal to 0.5.</exception>
		public VonMisesMaterial3D(double youngModulus, double poissonRatio, double yieldStress, double hardeningModulus)
		{
			this.youngModulus = youngModulus;

			if (poissonRatio == PoissonRatioForIncompressibleSolid)
			{
				throw new ArgumentException(
					"Poisson ratio cannot be" + PoissonRatioForIncompressibleSolid + "(incompressible solid)");
			}

			this.poissonRatio = poissonRatio;
			this.yieldStress = yieldStress;
			this.hardeningModulus = hardeningModulus;

			this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
			double lamda = (youngModulus * poissonRatio) / ((1 + poissonRatio) * (1 - (2 * poissonRatio)));
			double mi = youngModulus / (2 * (1 + poissonRatio));
			double value1 = (2 * mi) + lamda;

			this.elasticConstitutiveMatrix = Matrix.CreateZero(6, 6);
			this.elasticConstitutiveMatrix[0, 0] = value1;
			this.elasticConstitutiveMatrix[0, 1] = lamda;
			this.elasticConstitutiveMatrix[0, 2] = lamda;
			this.elasticConstitutiveMatrix[1, 0] = lamda;
			this.elasticConstitutiveMatrix[1, 1] = value1;
			this.elasticConstitutiveMatrix[1, 2] = lamda;
			this.elasticConstitutiveMatrix[2, 0] = lamda;
			this.elasticConstitutiveMatrix[2, 1] = lamda;
			this.elasticConstitutiveMatrix[2, 2] = value1;
			this.elasticConstitutiveMatrix[3, 3] = mi;
			this.elasticConstitutiveMatrix[4, 4] = mi;
			this.elasticConstitutiveMatrix[5, 5] = mi;

			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(PLASTIC_STRAIN, 0d),
				(STRESS_X, 0d),
				(STRESS_Y, 0d),
				(STRESS_Z, 0d),
				(STRESS_XY, 0d),
				(STRESS_XZ, 0d),
				(STRESS_YZ, 0d),
			});
		}

		public double[] Coordinates { get; set; }

		/// <summary>
		///   Gets the constitutive matrix.
		/// </summary>
		/// <value>
		///   The constitutive matrix.
		/// </value>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (this.constitutiveMatrix == null)
				{
					UpdateConstitutiveMatrixAndEvaluateResponse(new double[6]);
					this.constitutiveMatrix.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
				}
				return constitutiveMatrix;
			}
		}

		/// <summary>
		///   Gets the ID of the material.
		/// </summary>
		/// <value>
		///   The id.
		/// </value>
		public int ID => 1;

		/// <summary>
		///   Gets the incremental strains of the finite element's material.
		/// </summary>
		/// <value>
		///   The incremental strains.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
		/// </remarks>
		public double[] IncrementalStrains => this.incrementalStrains;

		/// <summary>
		///   Gets a value indicating whether this <see cref = "IStructuralMaterial" /> is modified.
		/// </summary>
		/// <value>
		///   <c>true</c> if modified; otherwise, <c>false</c>.
		/// </value>
		public bool IsCurrentStateDifferent() => modified;

		/// <summary>
		///   Gets the plastic strain.
		/// </summary>
		/// <value>
		///   The plastic strain.
		/// </value>
		public double PlasticStrain => this.plasticStrain;

		/// <summary>
		///   Gets the Poisson ratio.
		/// </summary>
		/// <value>
		///   The Poisson ratio.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
		/// </remarks>
		public double PoissonRatio
		{
			get
			{
				return this.poissonRatio;
			}
			set
			{
				this.poissonRatio = value;
			}
		}

		/// <summary>
		///   Gets the stresses of the finite element's material.
		/// </summary>
		/// <value>
		///   The stresses.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Stress_%28mechanics%29">Stress (mechanics)</a>
		/// </remarks>
		public double[] Stresses => this.stressesNew;

		/// <summary>
		///   Gets the Young's Modulus.
		/// </summary>
		/// <value>
		///   The young modulus.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
		/// </remarks>
		public double YoungModulus
		{
			get => this.youngModulus;
			set => this.youngModulus = value;
		}

		/// <summary>
		///   Creates a new object that is a copy of the current instance.
		/// </summary>
		/// <returns>
		///   A new object that is a copy of this instance.
		/// </returns>
		public object Clone()
		{
			return new VonMisesMaterial3D(this.youngModulus, this.poissonRatio, this.yieldStress, this.hardeningModulus)
			{
				modified = this.IsCurrentStateDifferent(),
				plasticStrain = this.plasticStrain,
				incrementalStrains = incrementalStrains.Copy(),
				stresses = stresses.Copy()
			};
		}

		/// <summary>
		///   Resets the indicator of whether the material is modified.
		/// </summary>
		public void ResetModified() => this.modified = false;

		/// <summary>
		///   Clears the stresses of the element's material.
		/// </summary>
		public void ClearStresses()
		{
			stresses.Clear();
			stressesNew.Clear();
		}

		public void ClearState()
		{
			modified = false;
			constitutiveMatrix.Clear();
			incrementalStrains.Clear();
			stresses.Clear();
			stressesNew.Clear();
			plasticStrain = 0;
			plasticStrainNew = 0;
		}

		/// <summary>
		///   Saves the state of the element's material.
		/// </summary>
		/// 
		public void SaveState()
		{
			this.plasticStrain = this.plasticStrainNew;
			stresses.CopyFrom(stressesNew);
		}
		public GenericConstitutiveLawState CreateState()
		{
			this.plasticStrain = this.plasticStrainNew;
			stresses.CopyFrom(stressesNew);
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(PLASTIC_STRAIN, plasticStrain),
				(STRESS_X, stresses[0]),
				(STRESS_Y, stresses[1]),
				(STRESS_Z, stresses[2]),
				(STRESS_XY, stresses[3]),
				(STRESS_XZ, stresses[4]),
				(STRESS_YZ, stresses[5]),
			});

			return currentState;
		}
		IHaveState ICreateState.CreateState() => CreateState();
		public GenericConstitutiveLawState CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				plasticStrain = currentState.StateValues[PLASTIC_STRAIN];
				stresses[0] = currentState.StateValues[STRESS_X];
				stresses[1] = currentState.StateValues[STRESS_Y];
				stresses[2] = currentState.StateValues[STRESS_Z];
				stresses[3] = currentState.StateValues[STRESS_XY];
				stresses[4] = currentState.StateValues[STRESS_XZ];
				stresses[5] = currentState.StateValues[STRESS_YZ];
			}
		}

		/// <summary>
		///   Updates the element's material with the provided incremental strains.
		/// </summary>
		/// <param name = "strainsIncrement">The incremental strains to use for the next step.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] strainsIncrement)
		{
			incrementalStrains.CopyFrom(strainsIncrement);
			this.CalculateNextStressStrainPoint();

			return stressesNew;
		}

		/// <summary>
		///   Builds the consistent tangential constitutive matrix.
		/// </summary>
		/// <param name = "value1"> This is a constant already calculated in the calling method. </param>
		/// <remarks>
		///   Refer to chapter 7.6.6 page 262 in Souza Neto. 
		/// </remarks>
		private Matrix BuildConsistentTangentialConstitutiveMatrix(double vonMisesStress, double[] unityvector)
		{
			Matrix consistenttangent = Matrix.CreateZero(TotalStresses, TotalStrains);
			this.constitutiveMatrix.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			consistenttangent.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			double dgamma = this.plasticStrainNew - this.plasticStrain;
			double v1 = -dgamma * 6 * Math.Pow(this.shearModulus, 2) / vonMisesStress;
			double Hk = 0;
			double Hi = this.hardeningModulus;
			double v2 = (dgamma / vonMisesStress - 1 / (3 * this.shearModulus + Hk + Hi)) * 6 * Math.Pow(this.shearModulus, 2);
			for (int i = 0; i < 6; i++)
			{
				for (int j = 0; j < 6; j++)
				{
					consistenttangent[i, j] = this.elasticConstitutiveMatrix[i, j] + v1 * SupportiveMatrixForConsistentConstitutiveMatrix[i, j] + v2 * unityvector[i] * unityvector[j];
				}
			}
			return consistenttangent;
		}

		/// <summary>
		///   Builds the tangential constitutive matrix.
		/// </summary>
		private void BuildTangentialConstitutiveMatrix()
		{
			this.constitutiveMatrix = Matrix.CreateZero(TotalStresses, TotalStrains);
			double invariantJ2New = this.GetDeviatorSecondStressInvariant(stressesNew);

			double value2 = (3 * this.shearModulus * this.shearModulus) /
							((this.hardeningModulus + (3 * this.shearModulus)) * invariantJ2New);

			var stressDeviator = this.GetStressDeviator(stressesNew);
			for (int k1 = 0; k1 < TotalStresses; k1++)
			{
				for (int k2 = 0; k2 < TotalStresses; k2++)
				{
					this.constitutiveMatrix[k2, k1] = this.elasticConstitutiveMatrix[k2, k1] -
													  (value2 * stressDeviator[k2] * stressDeviator[k1]);
				}
			}
		}

		/// <summary>
		///   Calculates the next stress-strain point.
		/// </summary>
		/// <exception cref = "InvalidOperationException"> When the new plastic strain is less than the previous one.</exception>
		private void CalculateNextStressStrainPoint()
		{
			var stressesElastic = new double[6];
			for (int i = 0; i < 6; i++)
			{
				stressesElastic[i] = this.stresses[i];
				for (int j = 0; j < 6; j++)
					stressesElastic[i] += this.elasticConstitutiveMatrix[i, j] * this.incrementalStrains[j];
			}

			double invariantJ2Elastic = this.GetDeviatorSecondStressInvariant(stressesElastic);
			double vonMisesStress = Math.Sqrt(3 * invariantJ2Elastic);
			double vonMisesStressMinusYieldStress = vonMisesStress -
													(this.yieldStress + (this.hardeningModulus * this.plasticStrain));

			bool materialIsInElasticRegion = vonMisesStressMinusYieldStress <= 0;

			if (materialIsInElasticRegion)
			{
				this.stressesNew = stressesElastic;
				this.constitutiveMatrix = this.elasticConstitutiveMatrix;
				this.plasticStrainNew = this.plasticStrain;
			}
			else
			{
				double deltaPlasticStrain = vonMisesStressMinusYieldStress /
											((3 * this.shearModulus) + this.hardeningModulus);
				this.plasticStrainNew = this.plasticStrain + deltaPlasticStrain;
				//double[] unityvector = GetStressDeviator(stressesNew);
				//invariantJ2Elastic = this.GetDeviatorSecondStressInvariant(stressesNew);
				double[] unityvector = GetStressDeviator(stressesElastic);
				invariantJ2Elastic = this.GetDeviatorSecondStressInvariant(stressesElastic);
				unityvector.ScaleIntoThis(Math.Sqrt(1 / (2 * invariantJ2Elastic)));
				for (int i = 0; i < 6; i++)
				{
					this.stressesNew[i] = stressesElastic[i] - (2 * this.shearModulus * deltaPlasticStrain * Math.Sqrt(1.5) * unityvector[i]);
				}
				this.constitutiveMatrix = this.BuildConsistentTangentialConstitutiveMatrix(vonMisesStress, unityvector);
				//this.BuildTangentialConstitutiveMatrix();
			}

			if (Math.Abs(this.plasticStrainNew) < Math.Abs(this.plasticStrain))
			{
				throw new InvalidOperationException("Plastic strain cannot decrease.");
			}

			this.modified = this.plasticStrainNew != this.plasticStrain;
		}

		/// <summary>
		///   Calculates and returns the first stress invariant (I1).
		/// </summary>
		/// <returns> The first stress invariant (I1).</returns>
		public double GetFirstStressInvariant(double[] stresses) => stresses[0] + stresses[1] + stresses[2];

		/// <summary>
		///   Calculates and returns the mean hydrostatic stress.
		/// </summary>
		/// <returns> The mean hydrostatic stress.</returns>
		public double GetMeanStress(double[] stresses) => GetFirstStressInvariant(stresses) / 3.0;

		/// <summary>
		///   Calculates and returns the second stress invariant (I2).
		/// </summary>
		/// <returns> The second stress invariant (I2).</returns>
		public double GetSecondStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1]) + (stresses[1] * stresses[2]) + (stresses[0] * stresses[2])
			- Math.Pow(stresses[5], 2) - Math.Pow(stresses[3], 2) - Math.Pow(stresses[4], 2);

		/// <summary>
		///   Calculates and returns the stress deviator tensor in vector form.
		/// </summary>
		/// <returns> The stress deviator tensor in vector form.</returns>
		public double[] GetStressDeviator(double[] stresses)
		{
			var hydrostaticStress = this.GetMeanStress(stresses);
			var stressDeviator = new double[]
			{
				stresses[0] - hydrostaticStress,
				stresses[1] - hydrostaticStress,
				stresses[2] - hydrostaticStress,
				stresses[3],
				stresses[4],
				stresses[5]
			};

			return stressDeviator;
		}

		/// <summary>
		///   Calculates and returns the third stress invariant (I3).
		/// </summary>
		/// <returns> The third stress invariant (I3). </returns>
		public double GetThirdStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1] * stresses[2]) + (2 * stresses[5] * stresses[3] * stresses[4])
			- (Math.Pow(stresses[5], 2) * stresses[2]) - (Math.Pow(stresses[3], 2) * stresses[0])
			- (Math.Pow(stresses[4], 2) * stresses[1]);

		/// <summary>
		///   Returns the first stress invariant of the stress deviator tensor (J1), which is zero.
		/// </summary>
		/// <returns> The first stress invariant of the stress deviator tensor (J1). </returns>
		public double GetDeviatorFirstStressInvariant(double[] stresses) => 0;

		/// <summary>
		///   Calculates and returns the second stress invariant of the stress deviator tensor (J2).
		/// </summary>
		/// <returns> The second stress invariant of the stress deviator tensor (J2). </returns>
		public double GetDeviatorSecondStressInvariant(double[] stresses)
		{
			double j2 = 0.0;
			j2 = 1.0 / 6.0 * (Math.Pow((stresses[0] - stresses[1]), 2) + Math.Pow((stresses[2] - stresses[1]), 2) + Math.Pow((stresses[0] - stresses[2]), 2));
			j2 += Math.Pow((stresses[3]), 2) + Math.Pow((stresses[4]), 2) + Math.Pow((stresses[5]), 2);
			return j2;
		}

		/// <summary>
		///   Calculates and returns the third stress invariant of the stress deviator tensor (J3).
		/// </summary>
		/// <returns> The third deviator stress invariant (J3). </returns>
		public double GetDeviatorThirdStressInvariant(double[] stresses)
		{
			double i1 = this.GetFirstStressInvariant(stresses);
			double i2 = this.GetSecondStressInvariant(stresses);
			double i3 = this.GetThirdStressInvariant(stresses);

			double j3 = (2 / 27 * Math.Pow(i1, 3)) - (1 / 3 * i1 * i2) + i3;
			return j3;
		}

	}
}
