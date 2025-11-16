package es.upm.dit.aled.lab4.genome;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Reads a FASTA file containing genetic information and allows for the search
 * of specific patterns within these data. The information is stored as an array
 * of bytes that contain nucleotides in the FASTA format. Since this array is
 * usually created before knowing how many characters in the origin FASTA file
 * are valid, an int indicating how many bytes of the array are valid is also
 * stored. All valid characters will be at the beginning of the array.
 * 
 * @author mmiguel, rgarciacarmona
 *
 */
public class FASTAReaderThreads {

	// All threads can access the same content and valid bytes since they are never
	// modified after the file is loaded.
	protected byte[] content;
	protected int validBytes;

	/**
	 * Creates a new FASTAReader from a FASTA file.
	 * 
	 * @param fileName The name of the FASTA file.
	 */
	public FASTAReaderThreads(String fileName) {
		try {
			this.readFile(fileName);
		} catch (IOException e) {
			System.out.println(e.getMessage());
			return;
		}
	}

	/*
	 * Helper method to read from a file. It populates the data array with upper
	 * case version of all the nucleotids found in the file. Throws an IOException
	 * if there is a problem accessing the file or the file is to big to fit in an
	 * array.
	 */
	private void readFile(String fileName) throws IOException {
		File f = new File(fileName);
		FileInputStream fis = new FileInputStream(f);
		DataInput fid = new DataInputStream(fis);
		long len = (int) fis.getChannel().size();
		if (len > Integer.MAX_VALUE) {
			fis.close();
			throw new IOException("The file " + fileName + " is too big. Can't be contained in an array.");
		}
		byte[] content = new byte[(int) len];
		int bytesRead = 0;
		int numRead = 0;
		String line;
		while ((line = fid.readLine()) != null) {
			// Put every character in upper case
			line = line.toUpperCase();
			numRead = line.length();
			byte[] newData = line.getBytes();
			for (int i = 0; i < numRead; i++)
				content[bytesRead + i] = newData[i];
			bytesRead += numRead;
		}
		fis.close();
		this.content = content;
		this.validBytes = bytesRead;
	}

	/**
	 * Provides the data array that contains the complete sequence of nucleotids
	 * extracted from the FASTA file.
	 * 
	 * @return The data array with each nucleotid taking one byte.
	 */
	public byte[] getContent() {
		return content;
	}

	/**
	 * Provides the amount of bytes in the data array that are valid. Since this
	 * array is created before the amount of bytes in the FASTA file that contain
	 * actual nucleotids are know, a worst-case scenario is assumed. So, only
	 * positions between content[0] and content[validBytes -1] have actual genomic
	 * data.
	 * 
	 * @return The number of valid bytes.
	 */
	public int getValidBytes() {
		return validBytes;
	}

	/**
	 * Realiza una búsqueda lineal para encontrar el patrón proporcionado en el array de datos.
	 * Para ello, crea un pool de hilos con tantos hilos como núcleos tenga el
	 * ordenador en el que se está ejecutando este código. Cada uno de estos hilos
	 * realizará una búsqueda lineal en el mismo byte[] (content), pero solo en un espacio
	 * de tamaño 1/(número de núcleos) del genoma, de modo que el trabajo se divide
	 * equitativamente entre todos los hilos.
	 * Cuando todos los hilos hayan terminado, agrega los resultados. Devuelve una lista de
	 * enteros que apuntan a las posiciones iniciales de todas las apariciones del
	 * patrón en los datos.
	 * 
	 * @param pattern El patrón que se desea encontrar.
	 * @return Todas las posiciones del primer carácter de cada aparición del
	 *         patrón en los datos.
	 */
	public List<Integer> search(byte[] pattern) { //busco ese pattern
		
		List<Integer> listaFinal = new ArrayList<Integer>(); //posiciones de dodne coincida la primera base del genoma para cada paattern que encuentre
		
		try {
			//Numero de cores que tiene mi ordenador
			int cores = Runtime.getRuntime().availableProcessors();
			System.out.println("Se usan " + cores + " cores");
			
			//Creo un executor --> jefe tareas: le mando las tareas a repartir, en vez de crear tantos threads como tareas
			ExecutorService executor = Executors.newFixedThreadPool(cores);
			
			/*esta listaFutures tiene "Future<List<Integer>>" --> cada cosa dentro de eso es un futuro;
			 * contiene una lista con las posiciones donde hay coincidencia--> para obtener cada lista
			 * es con: List<Integer> 
			 */
			List<Future<List<Integer>>> listaFutures = new ArrayList<>(); //fuera el for para guardarlos sin qu se reinicie en cada iteración
			
			int tamanoSegmento = content.length / cores;
			//Calculo los límites de cada segmento en el que trabaja cada core cada i en el for siguiente
			int lo = 0;
			int hi = tamanoSegmento;
			
			
			//cada i representa la tarea de cada core
			for(int i = 0; i < cores ; i++) {
				
				System.out.println("Creando tarea para core " + i + " con rango " + lo + " - " + hi);

				//Creo la tarea
				FASTASearchCallable tarea = new FASTASearchCallable(this, lo, hi, pattern);
				Future<List<Integer>> f = executor.submit(tarea); // aquí hay una tarea por core
				listaFutures.add(f);
				
				
				//Actualizo límites para el segmento siguiente
				lo += tamanoSegmento;
				hi += tamanoSegmento;
				
		}
			
			//Con este for, recorro listaFutures y hago que dejen de ser "futuros(promesas)" para que sean realidades --> .get()
			int indiceCore = 0; // contador fuera del for-each

			for (Future<List<Integer>> ff : listaFutures) {
			    List<Integer> real = ff.get();
			    System.out.println("La tarea del core " + indiceCore + " ha terminado, encontrando " + real.size() + " coincidencias");
			    listaFinal.addAll(real);
			    indiceCore++; // incrementa para el siguiente future
			}
			
			executor.shutdown();
			System.out.println("Búsqueda completa. Total coincidencias: " + listaFinal.size());

			
		} catch (Exception e) {
			System.out.println("Ha habido una interrupción " + e.getMessage());
		}
		
			
		return listaFinal;
		
	}

	public static void main(String[] args) {
		long t1 = System.nanoTime();
		FASTAReaderThreads reader = new FASTAReaderThreads(args[0]);
		if (args.length == 1)
			return;
		System.out.println("Tiempo de apertura de fichero: " + (System.nanoTime() - t1));
		long t2 = System.nanoTime();
		List<Integer> posiciones = reader.search(args[1].getBytes());
		System.out.println("Tiempo de búsqueda: " + (System.nanoTime() - t2));
		if (posiciones.size() > 0) {
			for (Integer pos : posiciones)
				System.out.println("Encontrado " + args[1] + " en " + pos);
		} else
			System.out.println("No he encontrado : " + args[1] + " en ningun sitio");
		System.out.println("Tiempo total: " + (System.nanoTime() - t1));
	}
}
