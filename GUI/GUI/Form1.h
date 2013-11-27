#pragma once

namespace GUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for Form1
	/// </summary>
	public ref class Form1 : public System::Windows::Forms::Form
	{
	public:
		Form1(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::TextBox^  in_unitcells_x;
	private: System::Windows::Forms::TextBox^  in_unitcells_y;
	private: System::Windows::Forms::TextBox^  in_unitcells_z;
	protected: 



	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::CheckBox^  in_pz;
	private: System::Windows::Forms::TextBox^  in_timestepsize;
	private: System::Windows::Forms::TextBox^  in_timesteps;
	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::ComboBox^  in_material;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::TextBox^  in_temp;
	private: System::Windows::Forms::Label^  label8;
	private: System::Windows::Forms::TextBox^  in_cutoff;
	private: System::Windows::Forms::Label^  label9;
	private: System::Windows::Forms::CheckBox^  in_thermostat;
	private: System::Windows::Forms::ProgressBar^  progressBar1;
	private: System::Windows::Forms::Button^  in_run;


	protected: 



	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->in_unitcells_x = (gcnew System::Windows::Forms::TextBox());
			this->in_unitcells_y = (gcnew System::Windows::Forms::TextBox());
			this->in_unitcells_z = (gcnew System::Windows::Forms::TextBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->in_pz = (gcnew System::Windows::Forms::CheckBox());
			this->in_timestepsize = (gcnew System::Windows::Forms::TextBox());
			this->in_timesteps = (gcnew System::Windows::Forms::TextBox());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->in_material = (gcnew System::Windows::Forms::ComboBox());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->in_temp = (gcnew System::Windows::Forms::TextBox());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->in_cutoff = (gcnew System::Windows::Forms::TextBox());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->in_thermostat = (gcnew System::Windows::Forms::CheckBox());
			this->progressBar1 = (gcnew System::Windows::Forms::ProgressBar());
			this->in_run = (gcnew System::Windows::Forms::Button());
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Times New Roman", 24, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label1->Location = System::Drawing::Point(141, 9);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(215, 36);
			this->label1->TabIndex = 0;
			this->label1->Text = L"MdSimulation";
			// 
			// in_unitcells_x
			// 
			this->in_unitcells_x->Location = System::Drawing::Point(12, 54);
			this->in_unitcells_x->Name = L"in_unitcells_x";
			this->in_unitcells_x->Size = System::Drawing::Size(100, 20);
			this->in_unitcells_x->TabIndex = 1;
			// 
			// in_unitcells_y
			// 
			this->in_unitcells_y->Location = System::Drawing::Point(12, 81);
			this->in_unitcells_y->Name = L"in_unitcells_y";
			this->in_unitcells_y->Size = System::Drawing::Size(100, 20);
			this->in_unitcells_y->TabIndex = 2;
			// 
			// in_unitcells_z
			// 
			this->in_unitcells_z->Location = System::Drawing::Point(13, 108);
			this->in_unitcells_z->Name = L"in_unitcells_z";
			this->in_unitcells_z->Size = System::Drawing::Size(100, 20);
			this->in_unitcells_z->TabIndex = 3;
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label2->Location = System::Drawing::Point(118, 55);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(94, 19);
			this->label2->TabIndex = 4;
			this->label2->Text = L"# Unitcells X";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label3->Location = System::Drawing::Point(118, 82);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(92, 19);
			this->label3->TabIndex = 5;
			this->label3->Text = L"# Unitcells Y";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label4->Location = System::Drawing::Point(118, 109);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(92, 19);
			this->label4->TabIndex = 6;
			this->label4->Text = L"# Unitcells Z";
			// 
			// in_pz
			// 
			this->in_pz->AutoSize = true;
			this->in_pz->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->in_pz->Location = System::Drawing::Point(13, 135);
			this->in_pz->Name = L"in_pz";
			this->in_pz->Size = System::Drawing::Size(185, 23);
			this->in_pz->TabIndex = 7;
			this->in_pz->Text = L"Periodic boundary in Z\?";
			this->in_pz->UseVisualStyleBackColor = true;
			// 
			// in_timestepsize
			// 
			this->in_timestepsize->Location = System::Drawing::Point(12, 165);
			this->in_timestepsize->Name = L"in_timestepsize";
			this->in_timestepsize->Size = System::Drawing::Size(100, 20);
			this->in_timestepsize->TabIndex = 8;
			// 
			// in_timesteps
			// 
			this->in_timesteps->Location = System::Drawing::Point(13, 192);
			this->in_timesteps->Name = L"in_timesteps";
			this->in_timesteps->Size = System::Drawing::Size(100, 20);
			this->in_timesteps->TabIndex = 9;
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label5->Location = System::Drawing::Point(118, 166);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(104, 19);
			this->label5->TabIndex = 10;
			this->label5->Text = L"Time step size";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label6->Location = System::Drawing::Point(118, 193);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(94, 19);
			this->label6->TabIndex = 11;
			this->label6->Text = L"# Time steps";
			// 
			// in_material
			// 
			this->in_material->DropDownStyle = System::Windows::Forms::ComboBoxStyle::DropDownList;
			this->in_material->FormattingEnabled = true;
			this->in_material->Location = System::Drawing::Point(13, 219);
			this->in_material->Name = L"in_material";
			this->in_material->Size = System::Drawing::Size(100, 21);
			this->in_material->TabIndex = 12;
			// Add things to the dropdownlist.
			this->in_material->Items->Add("Argon");

			
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label7->Location = System::Drawing::Point(118, 221);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(68, 19);
			this->label7->TabIndex = 13;
			this->label7->Text = L"Material";
			// 
			// in_temp
			// 
			this->in_temp->Location = System::Drawing::Point(12, 247);
			this->in_temp->Name = L"in_temp";
			this->in_temp->Size = System::Drawing::Size(100, 20);
			this->in_temp->TabIndex = 14;
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label8->Location = System::Drawing::Point(118, 248);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(95, 19);
			this->label8->TabIndex = 15;
			this->label8->Text = L"Temperature";
			// 
			// in_cutoff
			// 
			this->in_cutoff->Location = System::Drawing::Point(13, 274);
			this->in_cutoff->Name = L"in_cutoff";
			this->in_cutoff->Size = System::Drawing::Size(100, 20);
			this->in_cutoff->TabIndex = 16;
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label9->Location = System::Drawing::Point(118, 275);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(242, 19);
			this->label9->TabIndex = 17;
			this->label9->Text = L"Cuttoff multiples of lattice constant";
			// 
			// in_thermostat
			// 
			this->in_thermostat->AutoSize = true;
			this->in_thermostat->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->in_thermostat->Location = System::Drawing::Point(13, 300);
			this->in_thermostat->Name = L"in_thermostat";
			this->in_thermostat->Size = System::Drawing::Size(113, 23);
			this->in_thermostat->TabIndex = 18;
			this->in_thermostat->Text = L"Thermostat\?";
			this->in_thermostat->UseVisualStyleBackColor = true;
			// 
			// progressBar1
			// 
			this->progressBar1->Location = System::Drawing::Point(12, 426);
			this->progressBar1->Name = L"progressBar1";
			this->progressBar1->Size = System::Drawing::Size(494, 23);
			this->progressBar1->TabIndex = 19;
			// 
			// in_run
			// 
			this->in_run->Location = System::Drawing::Point(431, 374);
			this->in_run->Name = L"in_run";
			this->in_run->Size = System::Drawing::Size(75, 23);
			this->in_run->TabIndex = 20;
			this->in_run->Text = L"Run";
			this->in_run->UseVisualStyleBackColor = true;
			this->in_run->Click += gcnew System::EventHandler(this, &Form1::in_run_Click);
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(518, 461);
			this->Controls->Add(this->in_run);
			this->Controls->Add(this->progressBar1);
			this->Controls->Add(this->in_thermostat);
			this->Controls->Add(this->label9);
			this->Controls->Add(this->in_cutoff);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->in_temp);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->in_material);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->in_timesteps);
			this->Controls->Add(this->in_timestepsize);
			this->Controls->Add(this->in_pz);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->in_unitcells_z);
			this->Controls->Add(this->in_unitcells_y);
			this->Controls->Add(this->in_unitcells_x);
			this->Controls->Add(this->label1);
			this->Name = L"Form1";
			this->Text = L"Form1";
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	
	private: System::Void in_run_Click(System::Object^  sender, System::EventArgs^  e) {

				 // Start the simulation using all the parameters given in the gui.


			 }// end of Click
};
}

