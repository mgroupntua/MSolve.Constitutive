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
	public class VonMisesMaterial3DFunctionalHardening : IIsotropicContinuumMaterial3D
	{
		/// <summary>
		///   Identity second order tensor written in vector form.
		/// </summary>
		private static readonly double[] IdentityVector = new[] { 1.0, 1.0, 1.0, 0, 0, 0 };
		private const string EQUIVALENT_PLASTIC_STRAIN = "Equivalent plastic strain";
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
		///   Array for projecting a tensor on its deviatoric part.
		/// </summary>
		private static readonly double[,] DeviatoricProjection = new[,]
			{
				{  2.0/3.0, -1.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, 2.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, -1.0/3.0, 2.0/3.0, 0,   0,   0  },
				{  0,  0,  0, 1.0, 0,   0   },
				{  0,  0,  0, 0,   1.0, 0   },
				{  0,  0,  0, 0,   0,   1.0 }
			};

		/// <summary>
		///   Array for projecting a tensor on its volumetric part.
		/// </summary>
		private static readonly double[,] VolumetricProjection = new[,]
		{
			{ 1.0, 1.0,  1.0, 0,   0,   0   },
			{ 1.0, 1.0,  1.0, 0,   0,   0   },
			{ 1.0, 1.0,  1.0, 0,   0,   0   },
			{  0,  0,  0, 0, 0,   0   },
			{  0,  0,  0, 0,   0, 0   },
			{  0,  0,  0, 0,   0,   0}
		};

		/// <summary>
		///   The constitutive matrix of the material while still in the elastic region.
		/// </summary>
		private readonly Matrix elasticConstitutiveMatrix;

		/// <summary>
		///   Hardening modulus for linear hardening.
		/// </summary>
		private readonly Func<double, double> hardeningFunc;

		private double hardeningGrad;

		/// <summary>
		///   Hardening modulus for linear hardening.
		/// </summary>
		private double hardeningModulus;

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
		///   The bulk modulus.
		/// </summary>
		/// <remarks>
		///   <a href = "https://en.wikipedia.org/wiki/Bulk_modulus">Wikipedia - Bulk Modulus</a>
		/// </remarks>
		private readonly double bulkModulus;

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
		///   The current strain vector.
		/// </summary>
		private double[] strains = new double[6];

		/// <summary>
		///   The previously converged elastic strain vector.
		/// </summary>
		private double[] strainsElasticPrev;

		/// <summary>
		///   The current elastic strain vector.
		/// </summary>
		private double[] strainsElastic = new double[6];

		/// <summary>
		///   The previously converged plastic strain vector.
		/// </summary>
		private double[] strainsPlasticPrev;

		/// <summary>
		///   The current plastic strain vector.
		/// </summary>
		private double[] strainsPlastic = new double[6];

		/// <summary>
		///   The previously converged equivalent/accumulated plastic strain vector.
		/// </summary>
		private double strainsEquivalentPrev;

		/// <summary>
		///   The current equivalent/accumulated plastic strain vector.
		/// </summary>
		private double strainsEquivalent;

		/// <summary>
		///   The current increment stress vector.
		/// </summary>
		private double[] stresses = new double[6];

		/// <summary>
		///   The current iteration stress vector.
		/// </summary>
		private double[] stressesNew = new double[6];

		private bool matrices_not_initialized = true;

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
		/// <param name = "hardeningFunc">
		///   The hardening function.
		/// </param>
		/// <exception cref = "ArgumentException"> When Poisson ratio is equal to 0.5.</exception>
		public VonMisesMaterial3DFunctionalHardening(double youngModulus, double poissonRatio, double yieldStress, Func<double, double> hardeningFunc)
		{
			this.youngModulus = youngModulus;

			if (poissonRatio == PoissonRatioForIncompressibleSolid)
			{
				throw new ArgumentException(
					"Poisson ratio cannot be" + PoissonRatioForIncompressibleSolid + "(incompressible solid)");
			}

			this.poissonRatio = poissonRatio;
			this.yieldStress = yieldStress;
			this.hardeningFunc = hardeningFunc;

			this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
			this.bulkModulus = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double lamda = (youngModulus * poissonRatio) / ((1 + poissonRatio) * (1 - (2 * poissonRatio)));
			double mi = youngModulus / (2 * (1 + poissonRatio));
			double value1 = (2 * mi) + lamda;

			elasticConstitutiveMatrix = GetConstitutiveMatrix();
			InitializeMatrices();
		}

		public void InitializeMatrices()
		{
			strainsElasticPrev = new double[6];
			strainsPlasticPrev = new double[6];
			strainsEquivalentPrev = 0;
			matrices_not_initialized = false;
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
				if (this.constitutiveMatrix == null) UpdateConstitutiveMatrixAndEvaluateResponse(new double[6]);
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
		///   Gets a value indicating whether this <see cref = "IFiniteElementMaterial" /> is modified.
		/// </summary>
		/// <value>
		///   <c>true</c> if modified; otherwise, <c>false</c>.
		/// </value>
		public bool IsCurrentStateDifferent() => modified;

		/// <summary>
		///   Gets the plastic strain.
		/// </summary>
		/// <value>
		///   The plastic strain vector.
		/// </value>
		public double[] StrainsPlastic => this.strainsPlastic;

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
			return new VonMisesMaterial3DFunctionalHardening(this.youngModulus, this.poissonRatio, this.yieldStress, this.hardeningFunc)
			{
				modified = this.IsCurrentStateDifferent(),
				strainsEquivalent = this.strainsEquivalent,
				strainsEquivalentPrev = this.strainsEquivalentPrev,
				incrementalStrains = incrementalStrains.Copy(),
				stresses = stresses.Copy()
			};
		}

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
			incrementalStrains.Clear();
			stresses.Clear();
			stressesNew.Clear();
			strains.Clear();
			incrementalStrains.Clear();
			strainsElasticPrev.Clear();
			strainsElastic.Clear();
			strainsPlasticPrev.Clear();
			strainsPlastic.Clear();
			strainsEquivalent = 0;
			strainsEquivalentPrev = 0;
			constitutiveMatrix = GetConstitutiveMatrix();
		}

		/// <summary>
		/// Saves the state of the element's material.
		/// </summary>
		public GenericConstitutiveLawState CreateState()
		{
			stresses.CopyFrom(stressesNew);
			strainsElasticPrev.CopyFrom(strainsElastic);
			strainsPlasticPrev.CopyFrom(strainsPlastic);
			strainsEquivalentPrev = strainsEquivalent;
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(EQUIVALENT_PLASTIC_STRAIN, strainsEquivalent),
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
				strainsEquivalent = currentState.StateValues[EQUIVALENT_PLASTIC_STRAIN];
				stresses[0] = currentState.StateValues[STRESS_X];
				stresses[1] = currentState.StateValues[STRESS_Y];
				stresses[2] = currentState.StateValues[STRESS_Z];
				stresses[3] = currentState.StateValues[STRESS_XY];
				stresses[4] = currentState.StateValues[STRESS_XZ];
				stresses[5] = currentState.StateValues[STRESS_YZ];
			}
		}

		/// <summary>
		///   Resets the indicator of whether the material is modified.
		/// </summary>
		public void ResetModified() => this.modified = false;


		/// <summary>
		///   Updates the element's material with the provided incremental strains.
		/// </summary>
		/// <param name = "strainsIncrement">The incremental strains to use for the next step.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] strainsIncrement)
		{
			if (matrices_not_initialized)
			{ this.InitializeMatrices(); }
			incrementalStrains.CopyFrom(strainsIncrement);
			this.CalculateNextStressStrainPoint();

			currentState = new GenericConstitutiveLawState(this, new[]
{
				(EQUIVALENT_PLASTIC_STRAIN, strainsEquivalent),
				(STRESS_X, stresses[0]),
				(STRESS_Y, stresses[1]),
				(STRESS_Z, stresses[2]),
				(STRESS_XY, stresses[3]),
				(STRESS_XZ, stresses[4]),
				(STRESS_YZ, stresses[5]),
			});

			return stressesNew;
		}

		///// <summary>
		/////   Builds the consistent tangential constitutive matrix.
		///// </summary>
		///// <param name = "value1"> This is a constant already calculated in the calling method. </param>
		///// <remarks>
		/////   Refer to chapter 7.6.6 page 262 in Souza Neto. 
		///// </remarks>
		//private Matrix BuildConsistentTangentialConstitutiveMatrix(double vonMisesStress, double[] unityvector)
		//{
		//	Matrix consistenttangent = Matrix.CreateZero(TotalStresses, TotalStrains);
		//	double dgamma = this.strainsEquivalent - this.strainsEquivalentPrev;
		//	double v1 = -dgamma * 6 * Math.Pow(this.shearModulus, 2) / vonMisesStress;
		//	double Hk = 0;
		//	double Hi = this.hardeningRatio;
		//	double v2 = (dgamma / vonMisesStress - 1 / (3 * this.shearModulus + Hk + Hi)) * 6 * Math.Pow(this.shearModulus, 2);
		//	for (int i = 0; i < 6; i++)
		//	{
		//		for (int j = 0; j < 6; j++)
		//		{
		//			consistenttangent[i, j] = this.elasticConstitutiveMatrix[i, j] + v1 * DeviatoricProjection[i, j] + v2 * unityvector[i] * unityvector[j];
		//		}
		//	}
		//	return consistenttangent;
		//}

		/// <summary>
		///   Builds the consistent tangential constitutive matrix.
		/// </summary>
		/// <remarks>
		///   Refer to chapter 7.6.6 page 262 in Souza Neto. 
		/// </remarks>
		private void BuildConsistentTangentialConstitutiveMatrix(double[] strainDeviatoricTrial, double normStrainDeviatoricTrial, double vonMisesStress, double deltastrainsEquivalent)
		{
			var N = new double[incrementalStrains.Length];
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				N[i] = strainDeviatoricTrial[i] / normStrainDeviatoricTrial;
			}
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				for (int j = 0; j < incrementalStrains.Length; j++)
				{
					constitutiveMatrix[i, j] = 2 * shearModulus * (1 - (3 * shearModulus * vonMisesStress) / deltastrainsEquivalent) * DeviatoricProjection[i, j] +
						6 * shearModulus * shearModulus * (vonMisesStress / deltastrainsEquivalent - 1 / (3 * shearModulus + hardeningModulus)) * (N[i] * N[j]) + bulkModulus * VolumetricProjection[i, j];
				}
			}
		}

		///// <summary>
		/////   Builds the tangential constitutive matrix.
		///// </summary>
		//private void BuildTangentialConstitutiveMatrix()
		//{
		//	this.constitutiveMatrix = Matrix.CreateZero(TotalStresses, TotalStrains);
		//	double invariantJ2New = this.GetDeviatorSecondStressInvariant(stressesNew);

		//	double value2 = (3 * this.shearModulus * this.shearModulus) /
		//					((this.hardeningRatio + (3 * this.shearModulus)) * invariantJ2New);

		//	var stressDeviator = this.GetStressDeviator(stressesNew);
		//	for (int k1 = 0; k1 < TotalStresses; k1++)
		//	{
		//		for (int k2 = 0; k2 < TotalStresses; k2++)
		//		{
		//			this.constitutiveMatrix[k2, k1] = this.elasticConstitutiveMatrix[k2, k1] -
		//											  (value2 * stressDeviator[k2] * stressDeviator[k1]);
		//		}
		//	}
		//}

		/// <summary>
		///   Calculates the next stress-strain point.
		/// </summary>
		/// <exception cref = "InvalidOperationException"> When the new plastic strain is less than the previous one.</exception>
		private void CalculateNextStressStrainPoint()
		{
			//calculate trial variables
			var strainTrial = new double[incrementalStrains.Length];
			var stressTrial = new double[incrementalStrains.Length];
			var strainDeviatoricTrial = new double[incrementalStrains.Length];
			var stressVolumetricTrial = new double[incrementalStrains.Length];
			var stressDeviatoricTrial = new double[incrementalStrains.Length];
			var normStrainDeviatoricTrial = new double();
			for (int i = 0; i < incrementalStrains.Length; i++)
				strainTrial[i] = strainsElasticPrev[i] + incrementalStrains[i];
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				for (int j = 0; j < incrementalStrains.Length; j++)
				{
					strainDeviatoricTrial[i] = strainDeviatoricTrial[i] + DeviatoricProjection[i, j] * strainTrial[j];
					if (i >= 0 && i <= 2)
						stressTrial[i] += (2 * shearModulus * DeviatoricProjection[i, j] * strainTrial[j] + bulkModulus * VolumetricProjection[i, j] * strainTrial[j]);
					else
						stressTrial[i] += 0.5 * (2 * shearModulus * DeviatoricProjection[i, j] * strainTrial[j] + bulkModulus * VolumetricProjection[i, j] * strainTrial[j]);
				}
			}
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				for (int j = 0; j < incrementalStrains.Length; j++)
				{
					stressVolumetricTrial[i] += (double)1 / 3 * VolumetricProjection[i, j] * stressTrial[j];
					stressDeviatoricTrial[i] += DeviatoricProjection[i, j] * stressTrial[j];
				}
				normStrainDeviatoricTrial += strainDeviatoricTrial[i] * strainDeviatoricTrial[i];
			}
			normStrainDeviatoricTrial = Math.Sqrt(normStrainDeviatoricTrial);
			var J2 = GetDeviatorSecondStressInvariant(stressTrial);
			double vonMisesStress = Math.Sqrt(3 * J2);
			double vonMisesStressMinusYieldStress = vonMisesStress - hardeningFunc(strainsEquivalentPrev);

			bool materialIsInElasticRegion = vonMisesStressMinusYieldStress <= 0;
			var stressVolumetric = new double[incrementalStrains.Length];
			var stressDeviatoric = new double[incrementalStrains.Length];
			double deltastrainsEquivalent = 0;
			if (materialIsInElasticRegion)
			{
				stressesNew.CopyFrom(stressTrial);
				stressVolumetric.CopyFrom(stressVolumetricTrial);
				stressDeviatoric.CopyFrom(stressDeviatoricTrial);
				strainsEquivalent = strainsEquivalentPrev;
				constitutiveMatrix = elasticConstitutiveMatrix.CopyToFullMatrix();
			}
			else
			{
				strainsEquivalent = strainsEquivalentPrev;
				var yieldFunc = Math.Sqrt(3*J2) - hardeningFunc(strainsEquivalent);
				hardeningModulus = CalculateHardeningGradient(strainsEquivalent);
				while (Math.Abs(yieldFunc) >= 1e-6)
				{
					deltastrainsEquivalent += yieldFunc / (3 * shearModulus + hardeningModulus);
					strainsEquivalent = strainsEquivalentPrev + deltastrainsEquivalent;
					yieldFunc = Math.Sqrt(3 * J2) - 3 * shearModulus * deltastrainsEquivalent - hardeningFunc(strainsEquivalent);
					hardeningModulus = CalculateHardeningGradient(strainsEquivalent);
				}

				for (int i = 0; i < incrementalStrains.Length; i++)
				{
					stressVolumetric[i] = stressVolumetricTrial[i];
					stressDeviatoric[i] = (1 - (3 * shearModulus * deltastrainsEquivalent) / vonMisesStress) * stressDeviatoricTrial[i];
					stressesNew[i] = stressVolumetric[i] + stressDeviatoric[i];
				}
				BuildConsistentTangentialConstitutiveMatrix(strainDeviatoricTrial, normStrainDeviatoricTrial, deltastrainsEquivalent, vonMisesStress);
				//this.BuildTangentialConstitutiveMatrix();
			}

			if (Math.Abs(this.strainsEquivalent) < Math.Abs(this.strainsEquivalentPrev))
			{
				throw new InvalidOperationException("Plastic strain cannot decrease.");
			}

			this.modified = this.strainsEquivalent != this.strainsEquivalentPrev;

			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				if (i >= 0 && i <= 2)
					strainsElastic[i] = 1 / (2 * shearModulus) * stressDeviatoric[i] + 1 / (3 * bulkModulus) * stressVolumetric[i];
				else
					strainsElastic[i] = 2 * (1 / (2 * shearModulus) * stressDeviatoric[i] + 1 / (3 * bulkModulus) * stressVolumetric[i]);
				if (vonMisesStress > 0)
					strainsPlastic[i] = strainsPlasticPrev[i] + strainTrial[i] - strainsElastic[i];
				strains[i] = strainsElastic[i] + strainsPlastic[i];
			}
		}

		private double CalculateHardeningGradient(double strain_eq) => hardeningGrad = (hardeningFunc(strain_eq + 1e-8) - hardeningFunc(strain_eq)) / 1e-8;

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
			double i1 = this.GetFirstStressInvariant(stresses);
			double i2 = this.GetSecondStressInvariant(stresses);

			double j2 = (1 / 3d * Math.Pow(i1, 2)) - i2;
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

		private Matrix GetConstitutiveMatrix()
		{
			double fE1 = YoungModulus / (double)(1 + PoissonRatio);
			double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
			double fE3 = fE1 + fE2;
			double fE4 = fE1 * 0.5;
			var afE = Matrix.CreateZero(6, 6);
			afE[0, 0] = fE3;
			afE[0, 1] = fE2;
			afE[0, 2] = fE2;
			afE[1, 0] = fE2;
			afE[1, 1] = fE3;
			afE[1, 2] = fE2;
			afE[2, 0] = fE2;
			afE[2, 1] = fE2;
			afE[2, 2] = fE3;
			afE[3, 3] = fE4;
			afE[4, 4] = fE4;
			afE[5, 5] = fE4;

			return afE;
		}

	}
}
