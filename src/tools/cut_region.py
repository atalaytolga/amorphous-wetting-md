{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "59f16fa2",
   "metadata": {},
   "source": [
    "vapor simulation file: \"/home/mi/atalar93/data/simulation/vapor_thermalisation_T1_0_rho_0_1.h5\"\n",
    "\n",
    "liquid simulation file: \"/home/mi/atalar93/data/simulation/liquid_thermalisation_T1_0_rho_0_1.h5\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "8197ee8f",
   "metadata": {},
   "outputs": [],
   "source": [
    "import h5py\n",
    "import numpy as np"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "78498028",
   "metadata": {},
   "outputs": [],
   "source": [
    "# --- CONFIGURATION ---\n",
    "FILE_LIQ = \"/home/mi/atalar93/data/simulation/liquid_thermalisation_T0_7_rho_0_733.h5\"\n",
    "FILE_VAP = \"/home/mi/atalar93/data/simulation/vapor_thermalisation_T0_7_rho_0_0185.h5\"\n",
    "FILE_OUT = \"./out_hybrid_slab3.h5\"\n",
    "\n",
    "GAP_Z_WIDTH = 40.0\n",
    "BUFFER = 0.75  # Safety distance to prevent overlap at the interface"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "bccf85c5",
   "metadata": {},
   "outputs": [],
   "source": [
    "def wrap_pbc(positions, box_edges):\n",
    "    \"\"\"\n",
    "    Wraps coordinates into the centered box [-L/2, L/2].\n",
    "    \"\"\"\n",
    "    # Extract diagonal L from box edges\n",
    "    L = np.diag(box_edges)\n",
    "\n",
    "    return ((positions + L/2) % L) - L/2\n",
    "\n",
    "def remove_com_velocity(velocities):\n",
    "     return velocities - np.mean(velocities, axis = 0)\n",
    "\n",
    "\n",
    "def get_snapshot(file_handle):\n",
    "    \"\"\"\n",
    "    Extracts the last frame of position, velocity and the box from H5MD file.\n",
    "    \"\"\"\n",
    "\n",
    "    p_group = file_handle[\"particles\"][\"all\"] # Adjust group name if not 'all'\n",
    "\n",
    "    # Get the last frame [-1]\n",
    "\n",
    "    pos = p_group[\"position\"][\"value\"][-1]\n",
    "    vel = p_group[\"velocity\"][\"value\"][-1]\n",
    "    box = p_group[\"box\"][\"edges\"][()]\n",
    "\n",
    "    return pos, vel, box"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "edf870da",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Loading Data...\n"
     ]
    }
   ],
   "source": [
    "print(\"Loading Data...\")\n",
    "with h5py.File(FILE_LIQ, \"r\") as f_liq, h5py.File(FILE_VAP, 'r') as f_vap:\n",
    "    pos_liq, vel_liq, box_liq = get_snapshot(f_liq)\n",
    "    pos_vap, vel_vap, box_vap = get_snapshot(f_vap)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "19b4611a",
   "metadata": {},
   "outputs": [],
   "source": [
    "assert np.all(np.diag(box_vap)[:2] == np.diag(box_liq)[:2]), \"XY dimensions mismatch\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "0138447f",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([[ 40.,   0.,   0.],\n",
       "       [  0.,  40.,   0.],\n",
       "       [  0.,   0., 160.]])"
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "box_vap"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "41a6f9fd",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   Selected 23435 Liquid atoms\n",
      "   Selected 4092 Vapor atoms\n",
      "Writing ./out_hybrid_slab3.h5...\n",
      "Done\n"
     ]
    }
   ],
   "source": [
    "z_cut_min = - (GAP_Z_WIDTH / 2.0)\n",
    "z_cut_max = + (GAP_Z_WIDTH / 2.0)\n",
    "\n",
    "# Use Vapor box as Master\n",
    "L_master = np.diag(box_vap)\n",
    "\n",
    "pos_liq = wrap_pbc(pos_liq, box_vap)\n",
    "pos_vap = wrap_pbc(pos_vap, box_vap)\n",
    "\n",
    "z_vap_min = z_cut_min - BUFFER\n",
    "z_vap_max = z_cut_max + BUFFER\n",
    "\n",
    "mask_vap = (pos_vap[:, 2] < z_vap_min) | (pos_vap[:, 2] > z_vap_max)\n",
    "final_pos_vap = pos_vap[mask_vap]\n",
    "final_vel_vap = vel_vap[mask_vap]\n",
    "\n",
    "\n",
    "mask_liq = (pos_liq[:, 2] >= z_cut_min) & (pos_liq[:, 2] <= z_cut_max)\n",
    "final_pos_liq = pos_liq[mask_liq]\n",
    "final_vel_liq = vel_liq[mask_liq]\n",
    "\n",
    "\n",
    "print(f\"   Selected {len(final_pos_liq)} Liquid atoms\")\n",
    "print(f\"   Selected {len(final_pos_vap)} Vapor atoms\")\n",
    "\n",
    "\n",
    "# Concatenate Positions\n",
    "combined_pos = np.concatenate((final_pos_liq, final_pos_vap), axis=0)\n",
    "\n",
    "# Remove COM velocities:\n",
    "final_vel_liq = remove_com_velocity(final_vel_liq)\n",
    "final_vel_vap = remove_com_velocity(final_vel_vap)\n",
    "\n",
    "# Concatenate Velocities\n",
    "combined_vel = np.concatenate((final_vel_liq, final_vel_vap), axis=0)\n",
    "\n",
    "# Create Species Tags (Change np.zeros to np.ones if two species are required)\n",
    "species_liq = np.zeros(len(final_pos_liq), dtype=np.int32)\n",
    "species_vap = np.zeros(len(final_pos_vap), dtype=np.int32)\n",
    "combined_species = np.concatenate((species_liq, species_vap), axis=0)\n",
    "\n",
    "print(f\"Writing {FILE_OUT}...\")\n",
    "with h5py.File(FILE_OUT, 'w') as f_out:\n",
    "\n",
    "    # Copy H5MD Header from Liquid Simulation\n",
    "    with h5py.File(FILE_LIQ, \"r\") as source:\n",
    "        source.copy(\"h5md\", f_out)\n",
    "\n",
    "    #if \"parameters\" in source.keys():\n",
    "    #    source.copy(\"parameters\", f_out)\n",
    "\n",
    "    # Create particles group\n",
    "    p_all = f_out.create_group(\"particles/all\")\n",
    "\n",
    "    # Copy Positions (Time-dependent shape: [1, N, 3])\n",
    "    p_pos = p_all.create_group(\"position\")\n",
    "    p_pos.create_dataset(\"value\", data=combined_pos[None, ...], compression=\"gzip\")\n",
    "    p_pos[\"step\"] = np.array([0], dtype=np.int32)\n",
    "    p_pos[\"time\"] = np.array([0], dtype=np.int32)\n",
    "\n",
    "    # Copy Velocities (Time-dependent shape: [1, N, 3])\n",
    "    p_vel = p_all.create_group(\"velocity\")\n",
    "    p_vel.create_dataset(\"value\", data=combined_vel[None, ...], compression=\"gzip\")\n",
    "    p_vel[\"step\"] = np.array([0], dtype=np.int32)\n",
    "    p_vel[\"time\"] = np.array([0], dtype=np.int32)\n",
    "\n",
    "    # Create species \n",
    "    # TO-DO\n",
    "\n",
    "    # Create Box\n",
    "    p_box = p_all.create_group(\"box\")\n",
    "    p_box.create_dataset(\"edges\", data=box_vap)\n",
    "\n",
    "    with h5py.File(FILE_LIQ, \"r\") as src:\n",
    "        src_box = src[\"particles/all/box\"]\n",
    "        if \"boundary\" in src_box.attrs:\n",
    "            p_box.attrs[\"boundary\"] = src_box.attrs[\"boundary\"]\n",
    "        if \"dimension\" in src_box.attrs:\n",
    "            p_box.attrs[\"dimension\"] = src_box.attrs[\"dimension\"]\n",
    "print(\"Done\")\n",
    "    \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "7eeee8f2",
   "metadata": {},
   "outputs": [],
   "source": [
    "f_out.close()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "id": "74e52fcd",
   "metadata": {},
   "outputs": [],
   "source": [
    "vapor_width = 10\n",
    "cut_region = np.diag(box_liq).copy()\n",
    "cut_region[2] = vapor_width\n",
    "lower = -cut_region / 2.0\n",
    "upper = cut_region / 2.0"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 24,
   "id": "c2f99427",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[10. 10.  5.]\n",
      "[-10. -10.  -5.]\n"
     ]
    }
   ],
   "source": [
    "print(upper)\n",
    "print(lower)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "id": "ef6df1e9",
   "metadata": {},
   "outputs": [],
   "source": [
    "mask = np.all((pos_liq >= lower) & (pos_liq < upper), axis=1)\n",
    "idx = np.where(mask)[0]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c5846215",
   "metadata": {},
   "outputs": [],
   "source": [
    "def copy_particle_groups(file, inside):\n",
    "    shape = np.array([20, 20, 30])\n",
    "    lower = np.zeros_like(shape)\n",
    "    upper = lower + shape\n",
    "\n",
    "\n",
    "    padding = 0.75\n",
    "    # iterate over particle groups\n",
    "    for name, p in file[\"particles\"].items():\n",
    "        print(f\"group {name}:\", p)\n",
    "\n",
    "        # find particle indices that are within\n",
    "        L = np.diagonal(p[\"box/edges\"])\n",
    "        r = p[\"position/value\"][-1]\n",
    "\n",
    "        # Wrap positions inside the simulation cell\n",
    "        wrapped_pos = r - np.floor(r / L) * L\n",
    "        \n",
    "        mask = np.all((r >= lower) & (r < upper), axis=1)\n",
    "        if inside:\n",
    "            idx = np.where(mask)[0]\n",
    "        else:\n",
    "            idx = np.where(~mask)[0]\n",
    "\n",
    "        print (\"{0:d} particles selected\".format(idx.shape[0]))\n",
    "\n",
    "        pout = f_out.require_group(\"particles/\" + name)\n",
    "\n",
    "\n",
    "        for iname,item in p.items():\n",
    "            \n",
    "            # adjust box\n",
    "                if iname == \"box\":\n",
    "                    box = pout.require_group(\"box\")\n",
    "                    if \"edges\" not in box:\n",
    "                        box[\"edges\"] = np.diagflat(shape)\n",
    "\n",
    "                    boundary = item.attrs[\"boundary\"]\n",
    "                    boundary[np.where(shape != L)] = \"none\"\n",
    "                    box.attrs[\"boundary\"] = np.array(boundary.tolist(), dtype=h5py.special_dtype(vlen=bytes))\n",
    "                    box.attrs[\"dimension\"] = len(boundary)\n",
    "                    continue\n",
    "                group = pout.require_group(iname)\n",
    "                group[\"step\"] = [item[\"step\"][-1],]\n",
    "                group[\"time\"] = [item[\"time\"][-1],]\n",
    "                group.create_dataset(\n",
    "                    \"value\"\n",
    "                , data = [item[\"value\"][-1][idx],]\n",
    "                , compression=\"gzip\", compression_opts=6, shuffle=True\n",
    "                )"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "d183b008",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "group all: <HDF5 group \"/particles/all\" (3 members)>\n",
      "1872 particles selected\n",
      "group all: <HDF5 group \"/particles/all\" (3 members)>\n",
      "3797 particles selected\n"
     ]
    },
    {
     "ename": "TypeError",
     "evalue": "Can't broadcast (3797, 3) -> (1872, 3)",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mTypeError\u001b[0m                                 Traceback (most recent call last)",
      "Cell \u001b[0;32mIn [3], line 11\u001b[0m\n\u001b[1;32m      8\u001b[0m      f_liquid\u001b[38;5;241m.\u001b[39mcopy(\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mparameters\u001b[39m\u001b[38;5;124m\"\u001b[39m, f_out) \n\u001b[1;32m     10\u001b[0m copy_particle_groups(f_liquid, inside \u001b[38;5;241m=\u001b[39m \u001b[38;5;28;01mTrue\u001b[39;00m)\n\u001b[0;32m---> 11\u001b[0m \u001b[43mcopy_particle_groups\u001b[49m\u001b[43m(\u001b[49m\u001b[43mf_vapor\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43minside\u001b[49m\u001b[43m \u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43m \u001b[49m\u001b[38;5;28;43;01mFalse\u001b[39;49;00m\u001b[43m)\u001b[49m\n",
      "Cell \u001b[0;32mIn [2], line 66\u001b[0m, in \u001b[0;36mcopy_particle_groups\u001b[0;34m(file, inside)\u001b[0m\n\u001b[1;32m     64\u001b[0m append_ds(group, \u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mstep\u001b[39m\u001b[38;5;124m\"\u001b[39m, item[\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mstep\u001b[39m\u001b[38;5;124m\"\u001b[39m][\u001b[38;5;241m-\u001b[39m\u001b[38;5;241m1\u001b[39m])\n\u001b[1;32m     65\u001b[0m append_ds(group, \u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mtime\u001b[39m\u001b[38;5;124m\"\u001b[39m, item[\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mtime\u001b[39m\u001b[38;5;124m\"\u001b[39m][\u001b[38;5;241m-\u001b[39m\u001b[38;5;241m1\u001b[39m])\n\u001b[0;32m---> 66\u001b[0m \u001b[43mappend_ds\u001b[49m\u001b[43m(\u001b[49m\u001b[43mgroup\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[38;5;124;43mvalue\u001b[39;49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mitem\u001b[49m\u001b[43m[\u001b[49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[38;5;124;43mvalue\u001b[39;49m\u001b[38;5;124;43m\"\u001b[39;49m\u001b[43m]\u001b[49m\u001b[43m[\u001b[49m\u001b[38;5;241;43m-\u001b[39;49m\u001b[38;5;241;43m1\u001b[39;49m\u001b[43m]\u001b[49m\u001b[43m[\u001b[49m\u001b[43midx\u001b[49m\u001b[43m]\u001b[49m\u001b[43m)\u001b[49m\n",
      "Cell \u001b[0;32mIn [2], line 23\u001b[0m, in \u001b[0;36mcopy_particle_groups.<locals>.append_ds\u001b[0;34m(group, name, row)\u001b[0m\n\u001b[1;32m     21\u001b[0m dset \u001b[38;5;241m=\u001b[39m group[name]\n\u001b[1;32m     22\u001b[0m dset\u001b[38;5;241m.\u001b[39mresize((dset\u001b[38;5;241m.\u001b[39mshape[\u001b[38;5;241m0\u001b[39m] \u001b[38;5;241m+\u001b[39m \u001b[38;5;241m1\u001b[39m,) \u001b[38;5;241m+\u001b[39m dset\u001b[38;5;241m.\u001b[39mshape[\u001b[38;5;241m1\u001b[39m:])\n\u001b[0;32m---> 23\u001b[0m \u001b[43mdset\u001b[49m\u001b[43m[\u001b[49m\u001b[38;5;241;43m-\u001b[39;49m\u001b[38;5;241;43m1\u001b[39;49m\u001b[43m]\u001b[49m \u001b[38;5;241m=\u001b[39m row\n",
      "File \u001b[0;32mh5py/_debian_h5py_serial/_objects.pyx:54\u001b[0m, in \u001b[0;36mh5py._debian_h5py_serial._objects.with_phil.wrapper\u001b[0;34m()\u001b[0m\n",
      "File \u001b[0;32mh5py/_debian_h5py_serial/_objects.pyx:55\u001b[0m, in \u001b[0;36mh5py._debian_h5py_serial._objects.with_phil.wrapper\u001b[0;34m()\u001b[0m\n",
      "File \u001b[0;32m/usr/lib/python3/dist-packages/h5py/_debian_h5py_serial/_hl/dataset.py:980\u001b[0m, in \u001b[0;36mDataset.__setitem__\u001b[0;34m(self, args, val)\u001b[0m\n\u001b[1;32m    977\u001b[0m     mshape \u001b[38;5;241m=\u001b[39m val\u001b[38;5;241m.\u001b[39mshape\n\u001b[1;32m    979\u001b[0m \u001b[38;5;66;03m# Perform the write, with broadcasting\u001b[39;00m\n\u001b[0;32m--> 980\u001b[0m mspace \u001b[38;5;241m=\u001b[39m h5s\u001b[38;5;241m.\u001b[39mcreate_simple(\u001b[43mselection\u001b[49m\u001b[38;5;241;43m.\u001b[39;49m\u001b[43mexpand_shape\u001b[49m\u001b[43m(\u001b[49m\u001b[43mmshape\u001b[49m\u001b[43m)\u001b[49m)\n\u001b[1;32m    981\u001b[0m \u001b[38;5;28;01mfor\u001b[39;00m fspace \u001b[38;5;129;01min\u001b[39;00m selection\u001b[38;5;241m.\u001b[39mbroadcast(mshape):\n\u001b[1;32m    982\u001b[0m     \u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39mid\u001b[38;5;241m.\u001b[39mwrite(mspace, fspace, val, mtype, dxpl\u001b[38;5;241m=\u001b[39m\u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39m_dxpl)\n",
      "File \u001b[0;32m/usr/lib/python3/dist-packages/h5py/_debian_h5py_serial/_hl/selections.py:264\u001b[0m, in \u001b[0;36mSimpleSelection.expand_shape\u001b[0;34m(self, source_shape)\u001b[0m\n\u001b[1;32m    262\u001b[0m             eshape\u001b[38;5;241m.\u001b[39mappend(t)\n\u001b[1;32m    263\u001b[0m         \u001b[38;5;28;01melse\u001b[39;00m:\n\u001b[0;32m--> 264\u001b[0m             \u001b[38;5;28;01mraise\u001b[39;00m \u001b[38;5;167;01mTypeError\u001b[39;00m(\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mCan\u001b[39m\u001b[38;5;124m'\u001b[39m\u001b[38;5;124mt broadcast \u001b[39m\u001b[38;5;132;01m%s\u001b[39;00m\u001b[38;5;124m -> \u001b[39m\u001b[38;5;132;01m%s\u001b[39;00m\u001b[38;5;124m\"\u001b[39m \u001b[38;5;241m%\u001b[39m (source_shape, \u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39marray_shape))  \u001b[38;5;66;03m# array shape\u001b[39;00m\n\u001b[1;32m    266\u001b[0m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;28many\u001b[39m([n \u001b[38;5;241m>\u001b[39m \u001b[38;5;241m1\u001b[39m \u001b[38;5;28;01mfor\u001b[39;00m n \u001b[38;5;129;01min\u001b[39;00m remaining_src_dims]):\n\u001b[1;32m    267\u001b[0m     \u001b[38;5;66;03m# All dimensions from target_shape should either have been popped\u001b[39;00m\n\u001b[1;32m    268\u001b[0m     \u001b[38;5;66;03m# to match the selection shape, or be 1.\u001b[39;00m\n\u001b[1;32m    269\u001b[0m     \u001b[38;5;28;01mraise\u001b[39;00m \u001b[38;5;167;01mTypeError\u001b[39;00m(\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mCan\u001b[39m\u001b[38;5;124m'\u001b[39m\u001b[38;5;124mt broadcast \u001b[39m\u001b[38;5;132;01m%s\u001b[39;00m\u001b[38;5;124m -> \u001b[39m\u001b[38;5;132;01m%s\u001b[39;00m\u001b[38;5;124m\"\u001b[39m \u001b[38;5;241m%\u001b[39m (source_shape, \u001b[38;5;28mself\u001b[39m\u001b[38;5;241m.\u001b[39marray_shape))  \u001b[38;5;66;03m# array shape\u001b[39;00m\n",
      "\u001b[0;31mTypeError\u001b[0m: Can't broadcast (3797, 3) -> (1872, 3)"
     ]
    }
   ],
   "source": [
    "with h5py.File(\"/home/mi/atalar93/data/simulation/liquid_thermalisation_T1_0_rho_0_1.h5\", \"r\") as f_liquid, \\\n",
    "     h5py.File(\"/home/mi/atalar93/data/simulation/vapor_thermalisation_T1_0_rho_0_1.h5\", 'r') as f_vapor, \\\n",
    "     h5py.File(\"./out_test.h5\", 'w') as f_out:\n",
    "\n",
    "     # copy metadata\n",
    "     f_liquid.copy(\"h5md\", f_out)\n",
    "     if \"parameters\" in f_liquid.keys():\n",
    "          f_liquid.copy(\"parameters\", f_out) \n",
    "\n",
    "     copy_particle_groups(f_liquid, inside = True)\n",
    "     copy_particle_groups(f_vapor, inside = False)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "59a81bf5",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.2"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
