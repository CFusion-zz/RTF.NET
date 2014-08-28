using System;
using System.Diagnostics;
using System.IO;
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;

namespace RTF
{
	public class Window : GameWindow
	{
		public static Window active;
		private FField field;
		public Window()
			: base(512, 512,
			new GraphicsMode(), "RTF C#", 0,
			DisplayDevice.Default, 3, 0,
			GraphicsContextFlags.ForwardCompatible | GraphicsContextFlags.Debug)
		{
			active = this;
			field = new FField();
		}


		protected override void OnLoad(System.EventArgs e)
		{

			this.VSync = VSyncMode.Off;
		}

		protected override void OnUpdateFrame(FrameEventArgs e)
		{
			float dt = (float)((RenderTime + UpdateTime));
			this.Title = "FPS: " + 1 / dt;

			field.update(dt);

			if (Keyboard[OpenTK.Input.Key.Escape])
				Exit();
		}

		protected override void OnRenderFrame(FrameEventArgs e)
		{
			field.draw();
			SwapBuffers();
		}

		[STAThread]
		public static void Main()
		{
			using (Window wind = new Window())
			{
				wind.Run();
			}
		}
	}
}