using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.CompilerServices;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;

namespace RTF
{
	class FField
	{
		public int gridResolution = 128; // N
		public float dt = 0.1f; // deltaTime
		public float diffusionRate = 0.0f ;// diffusionRate;
		public float viscocity = 0.0f; // ~~
		public float force = 5.0f; // scales the mouse movement that generate a force
		public float source = 100.0f; // amount of density that will be deposited
		
		public float[] u;
 		public float[] v;
		public float[] dens;
		
		public float[] u_prev;
		public float[] v_prev;
		public float[] dens_prev;

		private int bufferSize;

		private int pmX = 0;
		private int pmY = 0;

		public FField()
		{
			bufferSize = (gridResolution + 2) * (gridResolution + 2); // Add an edge

			u = new float[bufferSize];
			v = new float[bufferSize];
			dens = new float[bufferSize];

			u_prev = new float[bufferSize];
			v_prev = new float[bufferSize];
			dens_prev = new float[bufferSize];

			Array.Clear(u, 0, bufferSize); // @32bit_warning
			Array.Clear(v, 0, bufferSize); // @32bit_warning
			Array.Clear(dens, 0, bufferSize); // @32bit_warning
			Array.Clear(u_prev, 0, bufferSize); // @32bit_warning
			Array.Clear(v_prev, 0, bufferSize); // @32bit_warning
			Array.Clear(dens_prev, 0, bufferSize); // @32bit_warning
		}
		private void handleInput()
		{
			int mX = Window.active.Mouse.X;
			int mY = Window.active.Mouse.Y;


			if (pmX == 0)
			{
				pmX = mX;
				pmY = mY;
			}

			if (!Window.active.Mouse[MouseButton.Left] && !Window.active.Mouse[MouseButton.Right])
				return;


			int i = (int)((mX / (float)512) * gridResolution + 1);
			int j = (int)(((512 - mY) / (float)512) * gridResolution + 1);

			if (i < 1 || i > gridResolution || j < 1 || j > gridResolution)
				return;

			if (Window.active.Mouse[MouseButton.Left])
			{
				u[IX(i, j)] = force * (mX - pmX);
				v[IX(i, j)] = force * (pmY - mY);
			}

			if (Window.active.Mouse[MouseButton.Right])
			{
				dens_prev[IX(i, j)] = source;
			}

			pmX = mX;
			pmY = mY;
		}

		public void update(float _dt)
		{
			dt = _dt;

			Array.Clear(u_prev, 0, bufferSize);
			Array.Clear(v_prev, 0, bufferSize);
			Array.Clear(dens_prev, 0, bufferSize);

			handleInput();
			
			velocityStep();
			densityStep();
		}

		public void draw()
		{
			GL.Viewport(0, 0, 512, 512);
			GL.MatrixMode(MatrixMode.Projection);
			GL.LoadIdentity();
			GL.Ortho(0.0, 1.0, 0.0, 1.0, -100.0, 100.0);
			GL.ClearColor(0.0f, 0.0f, 0.0f, 1.0f);
			GL.Clear(ClearBufferMask.ColorBufferBit);

			float h = 1.0f / gridResolution;

			GL.Begin(PrimitiveType.Quads);

			for (int i = 0; i <= gridResolution; i++)
			{
				float x = (i - 0.5f) * h;
				for (int j = 0; j <= gridResolution; j++)
				{
					float y = (j - 0.5f) * h;

					float d00 = dens[IX(i, j)];
					float d01 = dens[IX(i, j + 1)];
					float d10 = dens[IX(i + 1, j)];
					float d11 = dens[IX(i + 1, j + 1)];

					GL.Color3(d00, d00, d00); GL.Vertex2(x, y);
					GL.Color3(d10, d10, d10); GL.Vertex2(x + h, y);
					GL.Color3(d11, d11, d11); GL.Vertex2(x + h, y + h);
					GL.Color3(d01, d01, d01); GL.Vertex2(x, y + h);
				}
			}

			GL.End();

			/*GL.Color4(1.0f, 0.0f, 0.0f, 0.2f);
			GL.LineWidth(1.0f);

			GL.Begin(PrimitiveType.Lines);

			for (int i = 1; i <= field.gridResolution; i++)
			{
				float x = (i - 0.5f) * h;
				for (int j = 1; j <= field.gridResolution; j++)
				{
					float y = (j - 0.5f) * h;

					GL.Vertex2(x, y);
					GL.Vertex2(x + field.u[field.IX(i, j)] * 0.2, y + field.v[field.IX(i, j)] * 0.2);
				}
			}

			GL.End();
			*/
		}



		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public int IX(int i, int j)
		{
			return (i) + (gridResolution + 2) * (j);

		}

		void setBnd(int b, float[] x)
		{
			for (int i = 1; i <= gridResolution; i++)
			{
				x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
				x[IX(gridResolution + 1, i)] = b == 1 ? -x[IX(gridResolution, i)] : x[IX(gridResolution, i)];
				x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
				x[IX(i, gridResolution + 1)] = b == 2 ? -x[IX(i, gridResolution)] : x[IX(i, gridResolution)];
			}
			x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
			x[IX(0, gridResolution + 1)] = 0.5f * (x[IX(1, gridResolution + 1)] + x[IX(0, gridResolution)]);
			x[IX(gridResolution + 1, 0)] = 0.5f * (x[IX(gridResolution, 0)] + x[IX(gridResolution + 1, 1)]);
			x[IX(gridResolution + 1, gridResolution + 1)] = 0.5f * (x[IX(gridResolution, gridResolution + 1)] + x[IX(gridResolution + 1, gridResolution)]);
		}

		private void linearSolve(int b, float[] current, float[] prev, float a, float c)
		{
			for (int k = 0; k < 20; k++)
			{
				for (int i = 1; i <= gridResolution; i++)
				{
					for (int j = 1; j <= gridResolution; j++)
					{
						current[IX(i, j)] = (prev[IX(i, j)] + a *
							(current[IX(i - 1, j)] + current[IX(i + 1, j)] + current[IX(i, j - 1)] + current[IX(i, j + 1)]))
							/ c;
					}
				}
				setBnd(b, current);
			}
		}

		void advect(int b, float[] current, float[] prev, float[] u, float[] v)
		{
			int i0, j0, i1, j1;
			float x, y, s0, t0, s1, t1, dt0;

			dt0 = dt * gridResolution;
			for (int i = 1; i <= gridResolution; i++)
			{
				for (int j = 1; j <= gridResolution; j++)
				{
					x = i - dt0 * u[IX(i, j)];
					y = j - dt0 * v[IX(i, j)];

					if (x < 0.5f)
						x = 0.5f;

					if (x > gridResolution + 0.5f)
						x = gridResolution + 0.5f;
					i0 = (int)x; i1 = i0 + 1;

					if (y < 0.5f)
						y = 0.5f;

					if (y > gridResolution + 0.5f)
						y = gridResolution + 0.5f;

					j0 = (int)y;
					j1 = j0 + 1;
					s1 = x - i0;
					s0 = 1 - s1;
					t1 = y - j0;
					t0 = 1 - t1;

					current[IX(i, j)] = s0 * (t0 * prev[IX(i0, j0)] + t1 * prev[IX(i0, j1)]) +
								 s1 * (t0 * prev[IX(i1, j0)] + t1 * prev[IX(i1, j1)]);
				}
			}
			setBnd(b, current);
		}


		private void diffuse(int b, float[] current, float[] prev, float rate)
		{
			float a = dt * rate * gridResolution * gridResolution;
			linearSolve(b, current, prev, a, 1 + 4 * a);
		}



		private void project()
		{
			for (int i = 1; i <= gridResolution; i++)
			{
				for (int j = 1; j <= gridResolution; j++)
				{
					v_prev[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / gridResolution;
					u_prev[IX(i, j)] = 0;
				}
			}
			setBnd(0, v_prev);
			setBnd(0, u_prev);

			linearSolve(0, u_prev, v_prev, 1, 4);

			for (int i = 1; i <= gridResolution; i++)
			{
				for (int j = 1; j <= gridResolution; j++)
				{
					u[IX(i, j)] -= 0.5f * gridResolution * (u_prev[IX(i + 1, j)] - u_prev[IX(i - 1, j)]);
					v[IX(i, j)] -= 0.5f * gridResolution * (u_prev[IX(i, j + 1)] - u_prev[IX(i, j - 1)]);
				}
			}

			setBnd(1, u);
			setBnd(2, v);
		}

		private void densityStep()
		{
			addSource(dens, dens_prev);

			swapBuffers(ref dens_prev, ref dens);
			diffuse(0, dens, dens_prev, diffusionRate);

			swapBuffers(ref dens_prev, ref dens);
			advect(0, dens, dens_prev, u, v);
		}

		private void velocityStep()
		{
			// N, u, v, u_prev, v_prev, visc, dt
			addSource(u, u_prev);
			addSource(v, v_prev);

			swapBuffers(ref u_prev, ref u);
			diffuse(1, u, u_prev, viscocity);

			swapBuffers(ref v_prev, ref v);
			diffuse(2, v, v_prev, viscocity);

			project();

			swapBuffers(ref u_prev, ref u);
			swapBuffers(ref v_prev, ref v);

			advect(1, u, u_prev, u_prev, v_prev);
			advect(2, v, v_prev, u_prev, v_prev);

			project();
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private void swapBuffers(ref float[] a, ref float[] b)
		{
			float[] c;
			c = a;
			a = b;
			b = c;
		}

		private void addSource(float[] current, float[] prev)
		{
			// N X:u S:u_prev, dt
			for (int i = 0; i < bufferSize; i++)
			{
				current[i] += dt * prev[i];
			}
		}
	}
}
